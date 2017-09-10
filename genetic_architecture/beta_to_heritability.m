% Renaming of function 'beta_to_variance_explained'
% 
% Input: 
% beta - effect size 
% f - allele frequency
% V - total phenotypic variance 
% 
% Output: 
% h - variance explained 
% 
function h = beta_to_heritability(beta, f, V, SNP_TYPE, n_cases)

if(exist('n_cases', 'var')) % New: apply finite sample correction to get an unbiased estimator
    h = beta_to_variance_explained(beta, f, V, SNP_TYPE, n_cases);
else
    h = beta_to_variance_explained(beta, f, V, SNP_TYPE);
end
