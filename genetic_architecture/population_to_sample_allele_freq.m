% Take a random sample from a population and compute allele frequencies in sample
% Input:
% f_vec - frequency of each allele
% N - population size
% n_vec - number of individuals in sample for each allele (usually all the same)
%
% Output:
% k_vec - number of allele carriers for each allele (drawn at random)
%
function k_vec = population_to_sample_allele_freq(f_vec, N, n_vec, return_counts)

if(~exist('return_counts', 'var') || isempty(return_counts))
    return_counts = 1;
end

num_alleles = length(f_vec);
if(length(n_vec)==1)
    n_vec = repmat(n_vec, num_alleles, 1);
end
k_vec = zeros(size(n_vec));

for i=1:num_alleles % Sample without replacement
    k_vec(i) = hygernd(N, round(f_vec(i)*N), n_vec(i));
end
if(~return_counts) % return frequencies rather than counts (counts is defaule)
    k_vec = k_vec ./ n_vec;
end


