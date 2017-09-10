% Compute genotype matrix for a sample from a set of mutations and a coalescent tree
% Input: 
% pairwise_inds - indices representing a tree
% levels - coalescent level for each mutation 
% branches - branch of each mutation 
% 
% Output: 
% G - genotype matrix for the sample
% 
function G = genotype_mat_from_indices(pairwise_inds, levels, branches)

n = length(pairwise_inds)/2+1; 
n_sites = length(levels); 

G = zeros(n, n_sites); 

T = tree_from_pairwise_indices(pairwise_inds); 
for i=1:n_sites % get mutation carriers for each site 
    G(T{n+1-levels(i)}{branches(i)}, i) = 1; 
end


