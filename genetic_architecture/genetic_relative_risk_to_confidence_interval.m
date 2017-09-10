% Compute (approximate) st.d. for genetic-relative risk
% 
% Input: 
% grr_vec - genetic relative risk 
% f_vec - MAF
% mu_vec - prevalence 
% n_cases_vec - number of cases in GWAS study 
% n_controls_vec - number of controls in GWAS study 
% alpha - confidence level for declaring significance 
% 
% Output: 
% grr_std_vec - st.d. of estimated GRR 
% grr_min_vec - lower value of confidence interval for GRR 
% grr_max_vec - upper value of confidence interval for GRR 
% 
function [grr_std_vec grr_min_vec grr_max_vec] = genetic_relative_risk_to_confidence_interval( ...
    grr_vec, f_vec, mu_vec, n_cases_vec, n_controls_vec, alpha )

num_loci = length(f_vec);
p_z_x_marginal = ...
    genetic_relative_risk_to_p_z_x_marginal(f_vec, grr_vec, mu_vec);
p_z_x_marginal = pop_prob_to_case_control_prob(p_z_x_marginal, n_cases_vec, n_controls_vec); % Compute joint dist for case control 
%mu_vec = n_cases_vec ./ (n_cases_vec + n_controls_vec); % compute 'prevalence' in the stud

n_samples_vec = n_cases_vec + n_controls_vec; % get total number of samples 
theta = n_cases_vec ./ n_samples_vec; % get fraction of cases
contigency_table = p_z_x_marginal .* repmat(n_samples_vec, floor(num_loci / length(n_samples_vec)), 4);

old_var_estimator = 1; % so far first estimator looks better empirically 
if(old_var_estimator)
    grr_std_vec = sqrt(sum ( 1 ./ contigency_table, 2)); % compute standard deviation. The st.d. is on LOG scale!!!
else % new from Damjans Ph.D. thesis
    grr_std_vec = 1 ./ sqrt(2.* n_samples_vec .* f_vec .* (1-f_vec) .*  theta .* (1-theta));
end
    

if(~exist('alpha', 'var') || isempty(alpha)) % set confidence interval
    alpha = 0.05;
end
num_stds = norminv(1-alpha/2); % set # of st.d.s. to take

%grr_min_vec = grr_vec - num_stds.*grr_std_vec; % Compute 95% confidence interval
%grr_max_vec = grr_vec + num_stds.*grr_std_vec; % Compute 95% confidence interval
grr_min_vec = exp(log(grr_vec) - num_stds.*grr_std_vec); % Compute 95% confidence interval (use log)
grr_max_vec = exp(log(grr_vec) + num_stds.*grr_std_vec); % Compute 95% confidence interval (use log)
