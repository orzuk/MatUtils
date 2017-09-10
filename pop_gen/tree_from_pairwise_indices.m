

% Build tree from pairwise indices 
% Input: 
% pairwise_inds - indices of pairs of coalescing branches at each coalescent event
% 
% Output: 
% T - a tree represented as a cell array. Each entry T{i} corresponds to a level
% of the tree, and is a cell-array itself. Each sub-entry T{i}{j}
% represents a node and contains all descendent leafs of this nodes. 
% 
function T = tree_from_pairwise_indices(pairwise_inds)

n = length(pairwise_inds)/2+1; 

T = cell(n, 1); 
for i=1:n
    T{1}{i} = i; 
end
for j=1:n-1
   T{j+1} = T{j}; 
   min_ind = min(pairwise_inds((j*2-1):(j*2))); 
   max_ind = max(pairwise_inds((j*2-1):(j*2))); 
%   j_is = j
   T{j+1}{min_ind} = union(T{j+1}{pairwise_inds(j*2-1)}, T{j+1}{pairwise_inds(j*2)}); % add union
   T{j+1} = T{j+1}([1:(max_ind-1) (max_ind+1):end]); % remove 
end

