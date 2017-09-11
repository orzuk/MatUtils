% Compute sample size required to achieve a certain power.
% This is the 'inverse' of compute_association_power. Here we know the
% power and want to find the reuqired sample size
% Currently this is done by bisection method (we can optimize later)
%
% Input:
% p_mat - matrix of joint probabilites of genotype and phenotype (determines effect size)
% power_vec - vector with the power achieved
% n_cases_vec - vector of number of samples in test (cases or total).
%                   New:
%                   if this is in [0,1], it denote desired power. Then we need to compute the
%                   matching n!!!
% n_cases_frac - fraction of cases in a case-control study (default is 1/2).
% alpha_vec - probability of false posities (a set of p-value cutoffs)
% iters - number of simulations used to estimate p-value
% test_type - test performed (single locus, epsitasis, additive-model etc.)
% test_stat - statistical test used (hypergeometric, chi-square etc.)
% sampling_type - sample from the population, or do case-control balanced sampling (default, more realistic)
% const_effect_flag - this flag (default is 'OFF') says that when computing
%                     power, we simply assume that the input effect size was the observed one.
%                     Therefore, power is BINARY - either we passed the signficance threshold
%                     or not.
% model_params - additional parameters for test
% compute_sample_method - method for power computation: 'analytic' (default) or 'bisection'
%
% Output:
% n_samples_vec - vector of number of samples needed to achieve desired power.
% non_centrality_parameter - parameter capturing effect size (for analytic tests)
%
function [n_samples_vec non_centrality_parameter] = ...
    compute_sample_size_from_power(p_mat, power_vec, n_cases_frac, alpha_vec, ...
    iters, test_type, test_stat, sampling_type, const_effect_flag, model_params, compute_sample_method, varargin)


num_powers = length(power_vec); % number of different powers desired
if(~exist('n_cases_frac', 'var') || isempty(n_cases_frac))
    n_cases_frac = 0.5;
end
if(~exist('const_effect_flag', 'var'))
    const_effect_flag = [];
end
if(~exist('model_params', 'var'))
    model_params = [];
end
if(~exist('compute_sample_method', 'var'))
    compute_sample_method = 'analytic';
end
n_samples_vec = zeros(num_powers,1);

switch compute_sample_method % call compute power only once
    case 'analytic'
        if(~strcmp(test_stat, 'eric-crude-enrichment-analytic'))
            [~, ~, ~, non_centrality_parameter] = ...
                compute_association_power(p_mat, 1, 1, alpha_vec, ... % here make sure to adjust for case/controls
                iters, test_type, test_stat, sampling_type, const_effect_flag, model_params); % compute just the ncp
        end
end
for i=1:num_powers
    switch compute_sample_method
        case 'analytic' % assume chi-square distribution - always use the same test !!!
            switch strrep(test_stat, '_', '-')
                case 'eric-crude-enrichment-analytic'
                    [f_vec grr_vec mu_vec] = p_z_x_marginal_to_genetic_relative_risk(p_mat);
                    rr_vec = p_mat(:,4) ./ ( sum(p_mat(:, [2 4])) .* sum(p_mat(:, [3 4])) );                    
                    n_samples_vec(i) = 2 .* ( norminv(alpha_vec) - sqrt(rr_vec).*norminv(power_vec(i)) ) .^ 2 ./ ...
                        ((rr_vec-1).^2 .* f_vec);
                    
                otherwise
                    %            cur_NCP = fzero(@(x) ncx2cdf( chi2inv(1 - alpha_vec,1), 1, x) - 1 + power_vec(i), ...
                    %                chi2inv(1 - alpha_vec,1) * power_vec(i) * 2 );
                    cur_NCP = fminbnd(@(x) ( ncx2cdf( chi2inv(1 - alpha_vec,1), 1, x) - 1 + power_vec(i) ).^2, ...
                        0.1,  chi2inv(1 - alpha_vec,1) * 10 );
                    n_samples_vec(i) = round(2 .* cur_NCP / non_centrality_parameter(i));
            end % switch test stat
        case 'bisection'
            min_n_samples_vec = 2; max_n_samples_vec = 10000000;
            %    min_power = 0;     max_power = 1;
            mid_n_samples_vec = max_n_samples_vec/2;
            while ( (max_n_samples_vec - min_n_samples_vec) > 1)
                %    try_samples = mid_n_samples_vec
                n_cases_vec = mid_n_samples_vec * n_cases_frac;
                n_controls_vec = mid_n_samples_vec * (1-n_cases_frac);
                [mid_power, ~, ~, non_centrality_parameter] = ...
                    compute_association_power(p_mat, n_cases_vec, n_controls_vec, alpha_vec, ...
                    iters, test_type, test_stat, sampling_type, const_effect_flag, model_params);
                %    got_power = mid_power
                if(mid_power < power_vec(i)) % need to increase sample size
                    min_n_samples_vec = mid_n_samples_vec;
                    %            min_power = mid_power;
                else
                    max_n_samples_vec = mid_n_samples_vec;
                    %            max_power = mid_power;
                end
                mid_n_samples_vec = round((min_n_samples_vec+max_n_samples_vec)/2);
            end
            n_samples_vec(i) = mid_n_samples_vec;
    end % switch power method
end


