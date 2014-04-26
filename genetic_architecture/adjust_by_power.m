% Find power and add new loci
% 
% Input: 
% f_vec - MAF
% GRR_vec - genetic relative risk
% mu - prevalence
% n_cases - # of cases
% n_controls - # of controls 
%
% Output: 
% total_lambda_s - combined lambda_s of all loci (mult. model)
% total_lambda_s_adjusted - combined lambda_s of all loci adjusted for power
% total_h_liab - total heritability (liab. scale)
% total_h_liab_adjusted - total heritability (liab. scale) - power adjusted
% new_loci - number of predicted new loci 
%
function [total_lambda_s total_lambda_s_adjusted ...
    total_h_liab total_h_liab_adjusted new_loci] = ...
    adjust_by_power(f_vec, GRR_vec, mu, n_cases, n_controls)

alpha = 10^(-8); iters = 1000; % set power for detecting current loci 

p_z_x_marginal = genetic_relative_risk_to_p_z_x_marginal( f_vec, GRR_vec, mu);
loci_power_vec = compute_association_power(p_z_x_marginal, n_cases, n_controls, alpha, iters, ...
    'single-locus', 'chi-square', 'case-control'); % compute power to detect each locus
[lambda_s_vec lambda_s_add lambda_mz_add h_add V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab] = ...
    genetic_relative_risk_to_heritability(f_vec, GRR_vec, mu);

adjusted_lambda_s_vec = lambda_s_vec .^ (1./loci_power_vec); % adjust each lambda by power
total_lambda_s = prod(lambda_s_vec);
total_lambda_s_adjusted = prod(adjusted_lambda_s_vec);

