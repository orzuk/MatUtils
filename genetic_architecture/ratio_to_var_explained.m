% Compute variance explained by a genotype with a certain ratio
%
% Input: 
% R - ratio (enrichment)
% mu - prevalence 
% 
% Output: 
% V - variance explaiend (on disease scale) 
% 
function V = ratio_to_var_explained(R, mu)

f = 0.01; % population allele frequency 
V = genetic_relative_risk_to_variance_explained(f, R, mu, 'diploid');
