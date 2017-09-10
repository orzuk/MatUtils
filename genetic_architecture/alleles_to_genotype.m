% Collapse pairs of alleles to genotypes 
% 
% Input: 
% x_vec - vector of alleles (in adjacent pairs)
% x_probs - probability of each allele
% 
% Output: 
% g_vec - vector of genotypes (0,1,2 in each) 
% g_probs - probability of each genotype
% 
function [g_vec g_probs] = alleles_to_genotype(x_vec, x_probs)

g_vec = x_vec(:,1:2:end) + x_vec(:,2:2:end); % set genotype vector 
[g_vec , ~, J] = unique(g_vec, 'rows'); 
g_probs = accumarray(J, x_probs);
