% Take a random sample from a population and compute allele frequencies in sample
% Input:
% f_vec - frequency of each allele
% N - population size - total number of CHROMOSOMES !!!!
% n_vec - number of individual CHROMOSOMES!!! in sample for each allele (usually all the same)
% return_counts - flag saying to return counts (default) or frequencies 
%
% Output:
% k_vec - number of allele carriers for each allele (drawn at random)
%
function k_vec = population_to_sample_allele_freq(f_vec, N, n_vec, return_counts)

if(~exist('return_counts', 'var') || isempty(return_counts))
    return_counts = 1;
end
f_vec = vec2column(f_vec);
num_alleles = length(f_vec);
poiss_flag=0;
if(isscalar(n_vec))
    if(n_vec < 0.1*N)
        poiss_flag=1;
    end
    n_vec = repmat(n_vec, num_alleles, 1);
end

if((max(f_vec) <= 1.001) && ~isempty(find(f_vec(f_vec>0) < 1, 1))) % transfer alleles frequencies to counts
    f_vec = round(f_vec .* N);
end

if(poiss_flag) % might save time 
    poiss_inds = find(max(n_vec, round(f_vec)) < 0.1*N); hyg_inds = setdiff(1:num_alleles, poiss_inds);
    k_vec = zeros(num_alleles, 1); 
    k_vec(poiss_inds) = poissrnd( vec2column(round(f_vec(poiss_inds))) .* n_vec(poiss_inds) ./ N); % much faster: use poisson approximation
    k_vec(hyg_inds) = hygernd(N, vec2column(round(f_vec(hyg_inds))), n_vec(hyg_inds)); % Sample without replacement. Can sample together alleles with same frequencies? probably NOT!
else
    k_vec = hygernd(N, vec2column(round(f_vec)), n_vec); % Sample without replacement. Can sample together alleles with same frequencies? probably NOT!
end
if(~return_counts) % return frequencies rather than counts (counts is default)
    k_vec = k_vec ./ n_vec;
end


