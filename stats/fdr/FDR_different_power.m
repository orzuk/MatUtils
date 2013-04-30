% Compute the correct alphas for a tests with different statistical powers
% We start by doinga Bonferroni-like correction
% 
% Input: 
% n_vec - vector of sample sizes (influencing power)
% r - effect size (also influencing power)
% alpha - desried error probability (FWER) 
% 
% Output: 
% alpha_vec - vector of alpha threhsolds
% lambda - lagrange multiplier 
% optimal_power - total power (average number of rejections)
% naive_power - power when using the naive way (same alpha_i for all hypothesis)
%
function [alpha_vec lambda optimal_power naive_power] = FDR_different_power(n_vec, r, alpha) 

N = length(n_vec); 
lambda = fsolve(@(lambda) alpha_sum(lambda, n_vec, r, alpha), 0); 
alpha_vec = lambda_to_alpha(lambda, n_vec, r);

sum_sum_is = sum(alpha_vec)
check_solver_sum = alpha_sum(lambda, n_vec, r, alpha)

optimal_power = alpha_power(alpha_vec, n_vec, r)
naive_power = alpha_power(repmat(alpha/N, 1, N), n_vec, r)

function alpha_vec = lambda_to_alpha(lambda, n_vec, r) % We use lambda directly insread of log(lambda) for better numerics
alpha_vec = 1 - Phi( (2*lambda + r^2.*n_vec./2) ./ sqrt(2.*n_vec).*r );

function a_s = alpha_sum(lambda, n_vec, r, alpha)
a_s = sum(lambda_to_alpha(lambda, n_vec, r)) - alpha; 

function pow = alpha_power(alpha_vec, n_vec, r) % Compute statistical power (expected # of hypothesis rejected) for a given alpha vec

pow = sum(1 - Phi(  PhiInv(1-alpha_vec) - r .* sqrt(n_vec./2) ) );


