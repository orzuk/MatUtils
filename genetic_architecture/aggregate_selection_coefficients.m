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

rho_hat_s_population = corrcov(Sigma_hat_s_population); % compute corresponding esitmated correlation matrix 

% use estimated covariance to compute aggregated estimates
% Model is: for each gene g we have: 
% s_mat_g ~ N(s_g, s_population_Sigma_hat)

[num_genes, num_populations] = size(s_mat); 
s_combined = zeros(num_genes, 1); s_combined_std = zeros(num_genes, 1); 
for i=1:num_genes % loop on genes (each one has a slightly different covariance)
	
	Sigma_hat_cur_gene = s_mat_std(i,:) * rho_hat_s_population * s_mat_std(i,:)'; % first compute estimated covariance matrix for this gene 

	w = sum(inv(Sigma_hat_cur_gene)); w = w ./ sum(w); % compute weights vector 
	s_combined(i) = s_mat(i,:) * w; % combined estimator: linear combination  
	s_combined_std(i) = sqrt(w * Sigma_hat_cur_gene * w'); % combined estimator: estimated st.d. 	

end % loop on genes 
	

% Internal function for estimating standard deviation, from set of estimators for each gene 
function Sigma_hat_s_population = estimate_selection_cov_internal(s_mat, s_mat_std)

[num_genes, num_populations] = size(s_mat); 


% First estimate sigma 
s_vec_std = mean(s_mat_std); 
sigma_hat_vec = zeros(num_populations, 1); 
for j=1:num_populations
	sigma_hat_vec(j) = sqrt(var(s_mat(:,j)) - s_vec_std(j)^2); % take simple average of st.d. Should be replaced by MLE for each population!			
end
sigma_hat = mean(sigma_hat_vec); % estimator for width of TRUE selection coefficients distribution 

Sigma_hat_s_population = 0.5*diag(sigma_hat_vec.^2); % set variances (multiply by two later)
% Next estimate correlations
for j1=1:num_populations
	for j2=(j1+1):num_populations % loop on all pairs
		Sigma_hat_s_population(j1, j2) = 0.5 * ...
		(-var(s_mat(:,j1)-s_mat(:,j2)) + sigma_hat_vec(j1)^2+sigma_hat_vec(j2)^2);
	end
end
Sigma_hat_s_population = Sigma_hat_s_population+Sigma_hat_s_population'; % make symmetric
