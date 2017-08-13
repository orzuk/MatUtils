% Take a random sample from a population and compute allele frequencies in sample
% Input:
% x_vec - population frequency of alleles
% p_vec - population probability of being at this frequencies for alleles
% n - number of individuals in sample
% k_vec - (optional) evaluate sample SFS only for this set of values 
% cond_poly_flag - if 'on' we give conditional probabilities, given that allele is polymorphic
%
% Output:
% sample_x_vec - sample frequency of alleles
% sample_p_vec - sample probability of being at this frequencies for alleles
%
function [sample_x_vec, sample_p_vec] = population_to_sample_allele_freq_distribution( ...
    x_vec, p_vec, n, k_vec, cond_poly_flag) % n_vec, return_counts)

if(~exist('k_vec', 'var') || isempty(k_vec))
    k_vec = 0:n;
end
if(~exist('cond_poly_flag', 'var') || isempty(cond_poly_flag))
    cond_poly_flag = 0;
end

sample_x_vec = k_vec ./ n; % take possible allele frequencies in sample, including monomorphic alleles !
sample_p_vec = zeros(1, length(k_vec));

if(max(x_vec) >= 2) % transfer from counts to frequencies (must have highest frequency: 2N counts)
    x_vec = x_vec ./ max(x_vec); 
end

ctr=1;
for i= k_vec %  0:n % loop on allele frequencies
    % here use log-scale to avoid underflows
    sample_p_vec(ctr) = exp( log_binom(n, i) + my_log_sum_exp(i .* log(max(eps,x_vec)) + (n-i) .* log(max(eps,1-x_vec)) + log(p_vec)) ); ctr=ctr+1;
end

if(cond_poly_flag) % compute probability conditional on allelel polymorphic
    sample_p_vec = sample_p_vec ./ (sum(p_vec) - sum( exp(n .* log(max(eps,1-x_vec)) + log(p_vec)) ) - sum( exp(n .* log(max(eps,x_vec)) + log(p_vec)) ) );
    sample_p_vec(k_vec==0) = 0; sample_p_vec(k_vec == n) = 0; 
end
