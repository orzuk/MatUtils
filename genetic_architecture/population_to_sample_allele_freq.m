% Take a random sample from a population and compute allele frequencies in sample
% Input:
% f_vec - frequency of each allele
% N - population size - total number of CHROMOSOMES !!!! 
% n_vec - number of individual CHROMOSOMES!!! in sample for each allele (usually all the same)
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

if((max(f_vec) <= 1.001) && ~isempty(find(f_vec(f_vec>0) < 1, 1))) % transfer alleles frequencies to counts 
    f_vec = round(f_vec .* N); 
end

k_vec = hygernd(N, vec2column(round(f_vec)), n_vec); % Sample without replacement. Can sample together alleles with same frequencies? probably NOT!
% k_vec = zeros(size(n_vec));
% for i=1:num_alleles % Sample without replacement. Can sample together alleles with same frequencies? probably NOT!
%     k_vec(i) = hygernd(N, round(f_vec(i)), n_vec(i));  % slow part !!!! 
% end
if(~return_counts) % return frequencies rather than counts (counts is defaule)
    k_vec = k_vec ./ n_vec;
end


