% Compute the density of observed effect size, given a true effect size
%
% Input:
% true_beta  - true effect size (one value)
% f_vec - minor allele frequency 
% test_stat - both binary traits and QTL supported 
% alpha - significance level for declaring 'hit'
% beta_grid - values of effect size to consider
% num_cases - sample size (cases)
% num_controls - sample size (controls)
% study_type - discovery/replication/combined  
% simulation_flag - simulate or use gaussian approximation s
%
%
% Output:
% observed_beta_hist - histogram of observed effect size 
% observed_beta_mu - mean of observed effect size histogram
% observed_beta_std - standard deviation of observed effect size histogram
% 
function [observed_beta_hist observed_beta_mu observed_beta_std] = ...
    observed_effect_size_dist(true_beta, f_vec, test_stat, ...
    alpha, beta_grid, num_cases, num_controls, study_type, simulation_flag)


if(~exist('alpha', 'var') || isempty(alpha)) % set default significance level
    alpha = 5*10^(-8);
end
if(~exist('simulation_flag', 'var') || isempty(simulation_flag)) 
    simulation_flag = 0; % set default: analytic computation
end
if(~isempty(strfind(test_stat, 'QTL'))) % QTL
    qtl_flag = 1; trait_str = 'QTL';
    mu = [];
    observed_beta_mu = true_beta;
    observed_beta_std = 1 ./ sqrt(2 .* num_cases .* f_vec .*(1-f_vec)); % st.d. of beta
    test_type = 'chi-square-QTL';
else % binary trait
    qtl_flag = 0; trait_str = 'binary';
    true_grr = true_beta(1); mu = true_beta(2); true_beta = true_beta(1);  % encoding: first is GRR, than prevalence
    observed_beta_mu = true_beta(:,1);
    observed_beta_std = genetic_relative_risk_to_confidence_interval(true_grr, f_vec, mu, ...
        num_cases, num_controls); % compute st.d. for GRR estimation
    test_type = 'armitage';
end
    
iters = 10000; % set # of iterations for simulation 


switch simulation_flag % simulate loci 
    case 'simulate' % simulate data and take empirical distribution
        if(qtl_flag)
            p_vec = QTL_params_to_p_mat(f_vec, true_beta, 1);
        else
            p_vec = genetic_relative_risk_to_p_z_x_marginal(f_vec, true_grr, mu);
        end
        observed_p_vec = zeros(iters, length(p_vec));
        for i=1:iters
            if(mod(i, 100) == 0)
                simulate_iter = i
            end
            observed_p_vec(i,:) = ...
                simulate_genotype_phenotype(p_vec, ...
                num_cases, num_controls, 1, trait_str, 'allele'); % randomize beta according to 'true' randomization
        end
        
        if(qtl_flag) % get beta 
            observed_beta_vec = observed_p_vec(:,2);
        else % get GRR from case-control study
            [~, observed_beta_vec] = p_z_x_marginal_to_odds_ratio(...
                observed_p_vec ./ (num_cases+num_controls));
            CAF = observed_p_vec(:,3) ./ sum(observed_p_vec(:,[1 3]),2); % set control allele frequency  
            observed_beta_vec = odds_ratio_to_genetic_relative_risk(observed_beta_vec, ...
                CAF, mu); 
        end        
end

switch lower(study_type)
    case 'discovery' % discovery: assume censored gaussian distirbution
        threshold_beta = compute_discovery_boundary(f_vec, mu, num_cases, num_controls, ...
            alpha, 'marginal', test_type); % Determine discovery boundry
        switch simulation_flag
            case 'simulate'
                observed_beta_hist = ...
                    hist(observed_beta_vec(observed_beta_vec > threshold_beta), ...
                    beta_grid);
            otherwise % use gaussian approximation 
                observed_beta_hist = exp(-(beta_grid - true_beta).^2 ./ (2*observed_beta_std.^2)) ./ ...
                    (sqrt(2*pi) .* observed_beta_std); % symmetric gaussian
                observed_beta_hist(find(beta_grid < threshold_beta)) = 0; % truncate low effect sizes
                
        end
        
    case 'replication' % replication: assume gaussian distirbution
        switch simulation_flag
            case 'simulate' % simulate data and take empirical distribution
                observed_beta_hist = hist(observed_beta_vec, beta_grid);
            otherwise % use gaussian approximation 
                observed_beta_hist = exp(-(beta_grid - true_beta).^2 ./ (2*observed_beta_std.^2)) ./ ...
                    (sqrt(2*pi) .* observed_beta_std);
        end
    case 'combined'  % combined: mixture of discovery and replication
        observed_beta_hist = (num_cases(1) .* observed_effect_size_dist( true_beta, f_vec, test_stat, ...
            alpha, beta_grid, num_cases(1), num_controls, 'discovery', simulation_flag) + ...
            num_cases(2) .* observed_effect_size_dist( true_beta, f_vec, test_stat, ...
            alpha, beta_grid, num_cases(2), num_controls, 'replication', simulation_flag)) ./ ...
            (num_cases(1) + num_cases(2)); % take mixture (weighted average)
        return;
end
observed_beta_hist = normalize_hist(beta_grid, observed_beta_hist); % normalize distribution to sum to one
