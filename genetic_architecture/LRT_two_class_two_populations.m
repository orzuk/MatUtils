% Compute likelihood-ratio-test statistic between model 
% H0: both populations have same s value
% H1: populations have different s values 
% Statistic used: Likeihood-ratio-test statistic: 
% 
function LRT = LRT_two_class_two_populations(D1, k_vec1, n_vec1, LL1, D2, k_vec2, n_vec2, LL2)

% Assume that the parameters for a single population were already fitted 
LRT = 2*(LL1 + LL2); 


% Now fit together the data from both populations. Note: this is tough because for each one we've got 
% a different demography !!!! can we maximize together? depends how maximizing log-likelihood works. 
(Can be too slow to do grid-search) 
s_vec = logspace(???); 

LL_joints = ???; 


LRT = LRT - 2*LL_joint; 



% Aggregate alleles (should be transferred to a different function!)
% 
% Input: 
% s_mat - matrix of selection coefficients estimated for each population
% s_mat_std - matrix of standard deviations for selection coefficients estimated for each population
% 
% Output: 
% s_combined - combined estimator for selection
% s_combined_std - estimated standard deviation of combined estimator 
% 
function [s_combined, s_combined_std] = aggregate_selection_coefficients(s_mat, s_mat_std)

Sigma_hat_s_population = estimate_selection_cov_internal(s_mat, s_mat_std); % estimate covarinace matrix for s 

% use estimated covariance to compute aggregated estimates
% Model is: for each gene g we have: 
% s_mat_g ~ N(s_g, s_population_Sigma_hat)
s_combined = s_mat .* s_population_Sigma_hat 










