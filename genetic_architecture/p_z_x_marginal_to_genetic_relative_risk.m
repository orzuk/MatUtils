% Compute effect size and risk-allele-frequency
% from the joint distribution of genotype and phenotype
%
% Input:
% p_z_x_marginal - joint proability of genotype and phenotype (row vector of size 4 for each input)
%
% Output:
% f_vec - risk allele frequencies Pr(x=1)
% GRR_vec - genetic relative risk Pr(z=1|x=1) / Pr(z=1|x=0)
% mu_vec - prevalence Pr(z=1)
%
function [f_vec, GRR_vec, mu_vec] = p_z_x_marginal_to_genetic_relative_risk(p_z_x_marginal)

%n = size(p_z_x_marginal, 1); % number of different matrices 

f_vec = p_z_x_marginal(:,3) + p_z_x_marginal(:,4); % allele frequency
mu_vec = p_z_x_marginal(:,2) + p_z_x_marginal(:,4); % disease prevalence 

tmp_prob_one = p_z_x_marginal(:,4) ./ f_vec; % Pr(z=1 | x=1)
tmp_prob_zero = p_z_x_marginal(:,2) ./ (1-f_vec); % Pr(z=1 | x=0)
GRR_vec = tmp_prob_one ./ tmp_prob_zero; 

