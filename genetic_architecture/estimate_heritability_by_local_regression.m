% Estimate heritability by local linear regression near mean value of IBD 
%
% Input: 
% pairwise_IBD_mat - matrix of pairwise IBD levels for all individuals
% phenotype_vec - vector of phenotypes for each individual
% num_bins - number of bins for vizualization of histrogram
% plot_flag - 
% 
% Output: 
% h_all - estimated narrow sense heritability
% phenotype_corr_by_IBD_vec - structure containing phenotypc correlation statistics
% 
function [h_all phenotype_corr_by_IBD_vec] = estimate_heritability_by_local_regression(pairwise_IBD_mat, phenotype_vec, num_bins, plot_flag)

phenotype_vec = (phenotype_vec - mean(phenotype_vec)) ./ std(phenotype_vec); % Normalize
pairwise_IBD_mat = pairwise_IBD_mat - diag(diag(pairwise_IBD_mat)); % remove diagonal 

pairwise_IBD_vec = squareform(pairwise_IBD_mat); 
pairwise_phenotype_corr_mat = phenotype_vec * phenotype_vec';
pairwise_phenotype_corr_mat = pairwise_phenotype_corr_mat - diag(diag(pairwise_phenotype_corr_mat)); % remove diagonal
pairwise_phenotype_corr_vec = squareform(pairwise_phenotype_corr_mat); 


IBD_mean = mean(pairwise_IBD_mat(:)); 


IBD_estimation_range = [0, 2*IBD_mean]; % Take range centered at IBD_mean
f_range = find((pairwise_IBD_vec < IBD_estimation_range(2)) & ...
    (pairwise_IBD_vec >= IBD_estimation_range(1))); % get all IBD values within range

beta_linear_near_IBD = polyfit(pairwise_IBD_vec(f_range)', ...
    pairwise_phenotype_corr_vec(f_range)', 1); % local linear regression

h_all = beta_linear_near_IBD(1)*(1-IBD_mean)^2; % compute heritability estimator 

phenotype_corr_by_IBD_vec.IBD_bins = (0:num_bins) ./ num_bins;
phenotype_corr_by_IBD_vec.corr_mean = zeros(num_bins+1,1);
phenotype_corr_by_IBD_vec.corr_std = zeros(num_bins+1,1);
phenotype_corr_by_IBD_vec.corr_counts = zeros(num_bins+1,1);
cur_bins_vec = 1:max(floor(pairwise_IBD_vec*num_bins)+1);
phenotype_corr_by_IBD_vec.corr_mean(cur_bins_vec) = ...
    accumarray(floor(pairwise_IBD_vec'*num_bins)+1, pairwise_phenotype_corr_vec', [], @mean);
phenotype_corr_by_IBD_vec.corr_std(cur_bins_vec) = ...
    accumarray(floor(pairwise_IBD_vec'*num_bins)+1, pairwise_phenotype_corr_vec', [], @std);
phenotype_corr_by_IBD_vec.corr_counts(cur_bins_vec) = ...
    accumarray(floor(pairwise_IBD_vec'*num_bins)+1, pairwise_phenotype_corr_vec', [], @length);



