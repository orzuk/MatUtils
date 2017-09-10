% Combine different confidence intervals into one confidence interval 
% 
% Input: 
% mean_vec - vector of means for each separate estimator
% std_vec - vector of st.d. for each separate estimator
% corr_mat - matrix of pairwise correlations
% 
% Ouptut: 
% combined_mean - aggregate mean from all estimators 
% combined_std - aggregate st.d. from all estimator 
% 
function [combined_mean, combined_std] = CombineConfidenceIntervals(mean_vec, std_vec, corr_mat)


w = sum(inv(corr_mat)); w = w ./ sum(w); # get weights vector 
combined_mean = sum(w .* mean_vec); # get new combined estimator 

combined_std = sqrt( 1 ./ sum(sum(inv(corr_mat))) ); #  ???; % get st.d. of combined estimator ? need to check formula !!! 





% Estimate Covariance matrix from all populations: 
