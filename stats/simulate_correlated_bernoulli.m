% Simulates a set of correlated bernoulli random variables 
% with zero mean, unit variance and correlation coefficient rho.
% An alternative to sampleDichGauss01
%
% Input:
% n - number of vectors
% m - dimension
% rho - correlation coefficient
% f_vec - vector of means (optional. Default is zero)
% 
% Output
% x - random draws from a bernoulli 
%
function x = simulate_correlated_bernoulli(n, m, rho, f_vec)


x = simulate_correlated_gaussians(n, m, rho, zeros(m,1)); % simualte correlated latent gaussian r.v.s.
x_f_vec = norminv(1-f_vec); % get thresholds
x = x > repmat(x_f_vec, n, 1); % convert to bernoulli r.v.s. 

