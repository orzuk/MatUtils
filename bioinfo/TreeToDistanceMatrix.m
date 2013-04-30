% Compute the distance matrix associated with the tree 
% 
% Input: 
% Tree - a phylogenetic tree
% 
% The output: 
% D - a distance matrix between pairs of nodes
% 
function D = TreeToDistanceMatrix(Tree)

S = get(Tree, 'LEAFNAMES'); 
n = length(S); D = zeros(n); 

for i=1:n
    for j=i+1:n
        D(i,j) = GetBranchLength(Tree, {S{i}, S{j}}); 
    end
end

D = D+D'; 

