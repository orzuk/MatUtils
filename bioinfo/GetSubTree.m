% Extract a sub-tree from a subset of leaves of a phylogenetic tree
% 
% Input:
% Tree - the phylogenetic tree (can be input from file) 
% subspecies - vector of subspecies present in the tree
% 
% Output: 
% SubTree - a subtree containing only desired species
% 
function SubTree = GetSubTree(Tree, subspecies)

%if(isnumeric(subspecies))
%    subspecies = species_binary_to_vec(subspecies, Tree); 
%end
if(ischar(Tree)) % enable reading tree from file 
    Tree = phytreeread(Tree);
end
L = get(Tree, 'NUMLEAVES'); % get # of leaves
if(iscell(subspecies)) % here subspecies are given by names, and we transfer them to their indices
    LN = get(Tree, 'LEAFNAMES');
    [intersect_species, subspecies, J] = intersect(upper(LN), upper(subspecies)); % ignore case
end
if(length(subspecies) == L) % enable a binary 0/1 input
    if(max(subspecies) == 1)
        subspecies = find(subspecies);
    end
end
complement_subspecies = setdiff(1:L, subspecies);
SubTree = prune(Tree, complement_subspecies);

