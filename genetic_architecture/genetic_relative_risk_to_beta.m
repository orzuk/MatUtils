% Compute beta equivalent to a given GRR 
%
% Input:
% f_vec - risk allele frequencies
% GRR_marginal - genetic relative risk for each SNP
% mu - disease prevalence in population
% model - (optional) multiplicative or additive (or logistic?)
%
% Output:
% beta - effect size on liability scale
% beta_vec - individual effect sizes on liability scale
% 
function [beta  beta_vec] = ...
    genetic_relative_risk_to_beta(f_vec, GRR_marginal, mu)

[~,  h_liab_vec] = ...
    genetic_relative_risk_to_variance_explained(f_vec, GRR_marginal, mu, 'diploid');

beta_vec = heritability_to_beta(h_liab_vec, f_vec, 'diploid'); 
beta_vec = beta_vec .* ((-1).^(GRR_marginal < 1)); % take sign
beta = sum(beta_vec); % additive model 

