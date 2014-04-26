% Compute heritability from the joint distribution of genotype and phenotype
%
% Input:
% p_z_x_marginal - joint proability of genotype and phenotype (row vector of size 4 for each input)
% scale - (optional) scale of heritability. Default is liability 
% 
% Output:
% h_liab = heritability on the liability scale (default)
%
function h_liab = p_z_x_marginal_to_heritability(p_z_x_marginal, SNP_TYPE, scale)

if(~exist('SNP_TYPE', 'var') || isempty(SNP_TYPE))
    scale = 'diploid'; 
end
if(~exist('scale', 'var') || isempty(scale))
    scale = 'liability'; 
end

[f_vec, GRR_vec, mu_vec] = p_z_x_marginal_to_genetic_relative_risk(p_z_x_marginal);


h_liab = genetic_relative_risk_to_variance_explained(f_vec, GRR_vec, mu_vec, SNP_TYPE, scale);


