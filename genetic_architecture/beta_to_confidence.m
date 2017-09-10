% Compute (approximate) st.d. for beta parameter risk
function [beta_std_vec beta_min_vec beta_max_vec] = beta_to_confidence( ...
    beta_vec, f_vec,  n_samples_vec, alpha )

% Formula is: XXX, from http://en.wikipedia.org/wiki/Simple_linear_regression
beta_std_vec =  ( 1 ./ (f_vec.*n_samples_vec) +  1 ./ ((1-f_vec).*n_samples_vec)); % compute standard deviation.

if(~exist('alpha', 'var') || isempty(alpha)) % set confidence interval
    alpha = 0.05;
end
num_stds = norminv(1-alpha/2); % set # of st.d.s. to take
beta_min_vec = beta_vec - num_stds.*beta_std_vec; % Compute 95% confidence interval (use log)
beta_max_vec = beta_vec + num_stds.*beta_std_vec; % Compute 95% confidence interval (use log)
