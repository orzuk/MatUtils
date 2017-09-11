% Compute the odds ratio for phenotype from joint distribution of genotype and phenotype
%
% Input:
% p_z_x_marginal - joint proability of genotype and phenotype (row vector of size 4 for each input)
%
% Output:
% f_vec - risk allele frequencies Pr(x=1)
% OR_vec - genetic relative risk Pr(z=1|x=1) / Pr(z=1|x=0). This should be insensitive to sampling type!
% mu_vec - prevalence Pr(z=1)
%
function [f_vec, OR_vec, mu_vec] = p_z_x_marginal_to_odds_ratio(p_z_x_marginal)

%n = size(p_z_x_marginal, 1); % number of different matrices 

f_vec = p_z_x_marginal(:,3) + p_z_x_marginal(:,4); % allele frequency
mu_vec = p_z_x_marginal(:,2) + p_z_x_marginal(:,4); % disease prevalence 

OR_vec = p_z_x_marginal(:,1) .* p_z_x_marginal(:,4) ./ ...
    (p_z_x_marginal(:,2) .* p_z_x_marginal(:,3)); % multiply p(11)*p(00) / (p(01)*p(10))


