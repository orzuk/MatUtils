% Take a random sample from a population and compute allele frequencies in sample
function [k_vec, n_vec] = population_to_sample_allele_freq(f_vec, N, n_vec)

num_alleles = length(f_vec); 
if(length(n_vec)==1)
    n_vec = repmat(n_vec, num_alleles, 1); 
end
k_vec = zeros(size(n_vec)); 

% Sample without replacement 
for i=1:num_alleles
   k_vec(i) = hygernd(N, round(f_vec(i)*N), n_vec(i));      
end



