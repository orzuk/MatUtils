% Compute the total branch length of a subset of leaves of a phylogenetic tree
% 
% Input: 
% Tree - phylogeneteic tree (can be from file)
% subspecies - vector of species 
% 
% Output: 
% BL - branch length of the sub-tree containing species
%
function BL = GetBranchLength(Tree, subspecies, varargin)

if(exist('subspecies', 'var'))
    Tree = GetSubTree(Tree, subspecies);
end

B =  get(Tree, 'DISTANCES');
BL = sum(B); 
