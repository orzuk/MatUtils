% Take a random sample from a population and compute allele frequencies in sample
% Input:
% x_vec - population frequency of alleles
% p_vec - population probability of being at this frequencies for alleles
% n - number of individuals in sample
%
% Output:
% x_vec - sample frequency of alleles
% p_vec - sample probability of being at this frequencies for alleles
%
function [sample_x_vec, sample_p_vec] = population_to_sample_allele_freq_distribution(x_vec, p_vec, n) % n_vec, return_counts)

sample_x_vec = (0:n) ./ n; % take possible allele frequencies in sample, including monomorphic alleles !
sample_p_vec = zeros(1, n+1);


for i=0:n % loop on allele frequencies
    %    sample_p_vec(i+1) = nchoosek(n, i) * sum( x_vec.^i .* (1-x_vec).^(n-i) .* p_vec);
    
    % here use log-scale to avoid underflows
    sample_p_vec2(i+1) = exp( log_binom(n, i) + my_log_sum_exp(i .* log(max(eps,x_vec)) + (n-i) .* log(max(eps,1-x_vec)) + log(p_vec)) );
end



