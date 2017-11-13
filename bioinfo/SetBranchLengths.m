% Set the branch lengths of a subset of leaves of a phylogenetic tree
% 
% Input: 
% Tree - phylogeneteic tree (can be from file)
% BL - branch lengths 
% subspecies - vector of species 
% 
% Output: 
% NewTree - new tree with modified branch lengths
%
function NewTree = SetBranchLengths(Tree, BL)

D = get(Tree, 'Distances'); 
B = get(Tree, 'Pointers');
N = get(Tree, 'NodeNames');

NewTree = phytree(B, BL, N);

