% Compute similarity/distance between two phylogenetic trees 
% 
% Input: 
% tree1 - first tree
% w1 - first weights
% tree2 - second tree
% w2 - second weights
%
% The output: 
% S - a matrix of similarity measures 
% 
function S = TwoTreesToSimilarity(tree1, w1, tree2, w2)

bl1 = get(tree1, 'distances');
bl2 = get(tree2, 'distances'); % not used 

if(isempty(w1))
    w1 = get(tree1, 'distances'); 
end
if(isempty(w2))
    w2 = get(tree2, 'distances');
end
S = -sqrt(sum(bl1 .* (w1-w2).^2)); % For same topologies, take weighted euclidian distance


