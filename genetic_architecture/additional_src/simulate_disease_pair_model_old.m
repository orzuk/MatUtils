% Simulate disease probability according to any model with just two markers
%
% Input:
% f - genotypes vector
% rho_mat - matrix (2x2) giving disease probability for each of the four genotypes 
% d - disease contribution vector (additive)
% a - parameter of the disease prob. logistic function 1 / (1 +a*exp(-x))
% iters - number of simulations
%
% Output:
% rho - overall disease prob.
% V - variance in disease suseptability
% v_marginal - variance explained by each marker 
% v_both - variance explained by their combination 
%
function [rho V v_marginal v_both] = simulate_disease_pair_model_old(f, rho_mat)

% New: generic way (for any N)
N = length(f); 
bits_mat = zeros(2^N,N);
for i=1:N
    bits_mat(:,i) = bitget(0:2^N-1,i);
end
f_vec = prod(repmat(f,2^N,1) .^ (1-bits_mat),2) .* prod(repmat(1-f,2^N,1) .^ bits_mat,2); % generate a 2^N size vector

% First compute the frequency of alleles
%f_vec = [f(1)*f(2) f(1)*(1-f(2)) (1-f(1))*f(2) (1-f(1))*(1-f(2))]; 
%rho_vec = reshape(rho_mat, 1, 2^2);
rho = sum(f_vec .* rho_mat); % compute overall disease prob. 
V = rho*(1-rho); % compute overall variance 

rho_1 = zeros(N,1); rho_0 = zeros(N,1); 
for i=1:N
   rho_1(i) = sum(f_vec .* rho_mat .* (1-bits_mat(:,i))) ./ f(i);  
   rho_0(i) = sum(f_vec .* rho_mat .* bits_mat(:,i)) ./ (1-f(i));  
end

% rho_1_1 = f(2) * rho_mat(1,1) + (1-f(2)) * rho_mat(2,1);
% rho_1_0 = f(2) * rho_mat(1,2) + (1-f(2)) * rho_mat(2,2);
% rho_2_1 = f(1) * rho_mat(1,1) + (1-f(1)) * rho_mat(1,2);
% rho_2_0 = f(1) * rho_mat(2,1) + (1-f(1)) * rho_mat(2,2);
% v_marginal(1) = f(1) * rho_1_1 * (1-rho_1_1) + (1-f(1)) * rho_1_0 * (1-rho_1_0);
% v_marginal(2) = f(2) * rho_2_1 * (1-rho_2_1) + (1-f(2)) * rho_2_0 * (1-rho_2_0);

v_marginal = vec2column(f) .* rho_1 .* (1-rho_1) + (1-vec2column(f)) .* rho_0 .* (1-rho_0); 
v_both = sum(f_vec .* rho_mat .* (1-rho_mat));

