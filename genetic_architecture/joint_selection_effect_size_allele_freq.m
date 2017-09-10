% Compute joint and posterior probability of selection given allele frequency
%
% Input:
% s_bins - possible values of selection coefficient
% beta_bins - possible values of effect size 
% N - effective population size
% dist_params_struct - distribution-specific parameters structure
% coupling_str - type of joint distribution of effect size and selection coefficient
%
% Output:
% g_dist - joint distribution of effect size and selection coefficient
% s_hist - marginal distribution of selection coefficient
% beta_hist - marginal distribution of effect size 
%
function [g_hist s_hist beta_hist] = joint_selection_effect_size_allele_freq(s_bins, beta_bins, ... % allele frequency distribution
    N, dist_params_struct, coupling_str) % effective population size


num_s = length(s_bins); num_beta = length(beta_bins);
g_hist = zeros(num_s, num_beta);

% Compute weighted average of Hartl distribution, with weights given by the selection coefficient distirbution
switch coupling_str % allow different relationships between effect size and selection coefficient
    case 'eyre_walker' % use relationship from Eyre-Walker pnas 2010
        s_hist = gampdf(s_bins, dist_params_struct.k, dist_params_struct.theta); % First do the Gamma distribution
        mu = s_bins .^ dist_params_struct.tau;
        for i=1:num_s % Now compute the conditional distribution of beta given s
            g_hist(i,:) = s_hist(i) * 0.5 .* (...
                normpdf(beta_bins, mu(i), mu(i) * dist_params_struct.sigma) + ...
                normpdf(beta_bins, -mu(i), mu(i) * dist_params_struct.sigma)); % compute a Gaussian
        end
    case 'balancing' % assuming a balancing selection model
        beta_hist = normpdf(beta_bins, 0, 1); % assume standard Gaussian
        s_of_beta = exp(-beta_bins.^2 ./ (2.*(1+dist_params_struct.sigma^2))); % set deterministic relation
        for i=1:num_beta
            [~, cur_s_ind] = min(abs(s_bins - s_of_beta(i))); % Find closest index
            g_hist(cur_s_ind,i) = beta_hist(i);
        end
    case 'linear'

    case 'simple_eric' % simple case: 2 classes of mutations 
        s_hist(1) = 1 - dist_params_struct.alpha_null; % the proportion of mutations that are harmless
        [~, cur_s_ind] = min(abs(s_bins - dist_params_struct.null_s)); % Find closest index
        s_hist(cur_s_ind) = dist_params_struct.alpha_null; % null mutations 
        
        
        for i = [1 cur_s_ind] % no need to loop on everything 
            switch i
                case 1
                    [~, cur_beta_ind] = min(abs(beta_bins - dist_params_struct.null_beta)); % Find closest index
                otherwise
            end
            g_hist(i,cur_beta_ind) = s_hist(i); % set deterministic beta for each s 
        end % try different selection coefficients 
            
end % switch coupling type

g_hist = normalize_hist2d(s_bins, beta_bins, g_hist);
s_hist = sum(g_hist,2)'; beta_hist = sum(g_hist); % get marginal distributions 
s_hist = normalize_hist(s_bins, s_hist);
beta_hist = normalize_hist(beta_bins, beta_hist);
