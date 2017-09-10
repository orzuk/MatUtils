% Script for updating phenotypic correlation  counts for all iterations

if(~exist('all_iters_phenotype_corr_by_IBD_vec', 'var'))
    all_iters_phenotype_corr_by_IBD_vec.IBD_bins = phenotype_corr_by_IBD_vec.IBD_bins; 
    all_iters_phenotype_corr_by_IBD_vec.corr_mean = zeros(num_bins+1,1);
    all_iters_phenotype_corr_by_IBD_vec.corr_std = zeros(num_bins+1,1);
    all_iters_phenotype_corr_by_IBD_vec.corr_counts = zeros(num_bins+1,1);
end
all_iters_phenotype_corr_by_IBD_vec.corr_mean = all_iters_phenotype_corr_by_IBD_vec.corr_mean + ...
    phenotype_corr_by_IBD_vec.corr_mean .* phenotype_corr_by_IBD_vec.corr_counts; % Compute cumulative total
all_iters_phenotype_corr_by_IBD_vec.corr_std = all_iters_phenotype_corr_by_IBD_vec.corr_std + ...
    phenotype_corr_by_IBD_vec.corr_std.^2 .* phenotype_corr_by_IBD_vec.corr_counts + ...
    phenotype_corr_by_IBD_vec.corr_mean.^2 .* phenotype_corr_by_IBD_vec.corr_counts; % Compute cumulativevariance
all_iters_phenotype_corr_by_IBD_vec.corr_counts = all_iters_phenotype_corr_by_IBD_vec.corr_counts + ...
    phenotype_corr_by_IBD_vec.corr_counts; % Compute cumulative counts

