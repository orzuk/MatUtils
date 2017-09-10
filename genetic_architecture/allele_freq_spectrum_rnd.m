% Sample polymorphic alleles according to allele frequency distribution in a population at equilibrium. 
% Formula for density is from Sawyer&Hartl paper. We use 4*N*s (not 2*N*s) since N
% is the number of individuals (not chromosomes)
% The formula is: g(f) =  (1 - e^(s(1-q))) / [f(1-f)(1-e^(-S))]
%
% Input:
% s - selection coefficient (this is little s). Should be NEGATIVE for deleterious alleles (so fitness is 1+s)
% N - population size
% two_side_flag - return Derived (0) between 0 and 1 or Minor Allele Freq. (1) between 0 and 0.5
% M - num. alleles to sample
%
% Output:
% f - allele frequencies of sampled alleles
%
function f = allele_freq_spectrum_rnd(s, N, two_side_flag, M)

x_vec = (1:2*N-1) ./ (2*N); p_vec = exp(allele_freq_spectrum(x_vec, s, N, two_side_flag, 'log')); % compute density

f = weighted_rand(p_vec, M) ./ (2.*N);

