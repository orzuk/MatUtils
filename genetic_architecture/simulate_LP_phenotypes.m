% Simulate phenotype according to LP model
function LP_phenotype_vec = simulate_LP_phenotypes(SNP_mat, f_vec, num_pathways, h_pathway, h_shared_env, trait_type, mu)

[num_snps, num_individuals] = size(SNP_mat);

[snp_start_inds snp_end_inds] = ...
    divide_region_to_blocks(1, num_snps, ceil(num_snps/num_pathways));
beta_snps = sqrt(num_pathways / sum(f_vec .* (1-f_vec))); % set coefficients such that variance sums to one

SNP_mat = SNP_mat - repmat(f_vec, 1, num_individuals); % normalize genotypes to have mean zero


LP_phenotype_vec = zeros(num_individuals,1)-99999999999; % start with very negative
for i=1:num_pathways
    temp_pathay_vec = ...
        sum(beta_snps .* SNP_mat(snp_start_inds(i):snp_end_inds(i),:))' .* sqrt(h_pathway) + ...
        randn(num_individuals,1) .* sqrt(1-h_pathway);% Generate phenotype
    LP_phenotype_vec = max(LP_phenotype_vec, temp_pathay_vec);
end

switch trait_type
    case {'disease', 'binary'}
        x_mu = norminv(1-mu); % set threshold for disease
        LP_phenotype_vec = LP_phenotype_vec > x_mu;
end
