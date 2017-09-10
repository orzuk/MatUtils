% Compute the joint distribution of genotype and phenotype from effect size
% and risk-allele-frequency
%
% Input:
% f_vec - risk allele frequencies Pr(x=1) IN THE POPULATION (chenges if we move to case-control)
% GRR_vec - genetic relative risk Pr(z=1|x=1) / Pr(z=1|x=0)
% mu_vec - prevalence Pr(z=1)
%
% Output:
% p_z_x_marginal - joint proability of genotype and phenotype (row vector of size 4 for each input)
%
function p_z_x_marginal = genetic_relative_risk_to_p_z_x_marginal(f_vec, GRR_vec, mu_vec)

p_z_x_marginal = zeros(max(length(f_vec), length(GRR_vec)), 4);
f_vec = vec2column(f_vec);
GRR_vec = vec2column(GRR_vec);
mu_vec = vec2column(mu_vec);

tmp_prob_zero = mu_vec ./ (1-f_vec + GRR_vec .* f_vec);
tmp_prob_one = GRR_vec .* tmp_prob_zero;

p_z_x_marginal(:,1) = (1-f_vec) .* (1-tmp_prob_zero); % Pr(x_i=0,z=0)
p_z_x_marginal(:,2) = (1-f_vec) .* tmp_prob_zero; % Pr(x_i=0,z=1)
p_z_x_marginal(:,3) = f_vec .* (1-tmp_prob_one); % Pr(x_i=1,z=0)
p_z_x_marginal(:,4) = f_vec .* tmp_prob_one; % Pr(x_i=1,z=1)

