% Extract a sub-tree from a subset of leaves of a phylogenetic tree
% 
% Input; 
% Tree - the phylogenetic tree (can be input from file) 
% output_dir - where to save all tree files 
% 
% Output: 
% SubTrees - a structure with all possible subtrees obtained by cutting the tree at one point in a branch 
% 
function AllSubTrees = GetAllConnectedSubTrees(Tree, output_dir)

%if(isnumeric(subspecies))
%    subspecies = species_binary_to_vec(subspecies, Tree); 
%end
if(ischar(Tree)) % enable reading tree from file 
    tree_header = loadcellfile(Tree); % read also .text file 
    tree_line = strmatch('TREE:', tree_header(:,1)); 

    Tree = phytreeread(Tree);
end


num_leaves = get(Tree, 'numleaves'); % get # of leaves
leaf_names = get(Tree, 'leafnames');
node_names = get(Tree, 'nodenames'); 
num_nodes = get(Tree, 'numnodes'); 

sub_trees = cell(num_nodes+num_leaves); 
for i=num_leaves+1:num_nodes % loop only on non-leaves 
    sub_trees{i} = subtree(Tree,i); % node_names{i})
end

for i=num_leaves+1:num_nodes % loop only on non-leaves % Get complement trees of all trees (list them first)
    subtree_leaf_names = get(sub_trees{i}, 'nodenames');
    complement_subtree_leaf_names = setdiff(leaf_names, subtree_leaf_names);
    if(length(complement_subtree_leaf_names) > 1) % take only trees with at least 2 nodes
        sub_trees{i-num_leaves} = GetSubTree(Tree, complement_subtree_leaf_names);
    end
end    
for i=num_nodes+1:(num_nodes+num_leaves) % add all complements of one species 
    complement_subtree_node_names = setdiff(node_names, node_names{i-num_nodes});
    sub_trees{i} = GetSubTree(Tree, complement_subtree_node_names); 
end
    
keep_inds_vec = zeros(num_nodes+num_leaves,1); 
node_names_vec = cell(num_nodes+num_leaves,1);
leaf_names_vec = cell(num_nodes+num_leaves,1); 
binary_vec = zeros(num_nodes+num_leaves, 1, 'uint64');
for i=1:num_nodes+num_leaves % remove empty trees and trees with one node
    if(isempty(sub_trees{i}))
        continue;
    end
    node_names_vec{i} = get(sub_trees{i}, 'nodenames');
    leaf_names_vec{i} = get(sub_trees{i}, 'leafnames'); 
    if(get(sub_trees{i}, 'numnodes') > 1)
        keep_inds_vec(i) = 1; 
    end        
    [~, I] = intersect(leaf_names, leaf_names_vec{i});  % Get binary representation
    for j=I' % loop on species 
        binary_vec(i) = bitset(binary_vec(i), j);
    end
end
    
sub_trees = sub_trees(find(keep_inds_vec)); 
%node_names_vec = node_names_vec(find(keep_inds_vec)); 
binary_vec = binary_vec(find(keep_inds_vec)); 

[~, keep_inds] = unique(binary_vec);
AllSubTrees = sub_trees(keep_inds); 
binary_vec = binary_vec(keep_inds); 

if(exist('output_dir', 'var') && (~isempty(output_dir))) % save to files 
    for i=1:length(AllSubTrees)
        my_mkdir(fullfile(output_dir, ['clade_' num2str(binary_vec(i))])); 
        cur_output_file = fullfile(output_dir, ...
            ['clade_' num2str(binary_vec(i))], ['tree_clade_' num2str(binary_vec(i)) '.mod']);
        phytreewrite(cur_output_file, AllSubTrees{i}, 'branchnames',false); 
        tree_cell = loadcellfile(cur_output_file); tree_cell{1,1} = ['TREE: ' tree_cell{1,1}]; 
        tree_cell_combined = tree_header(1:tree_line-1,:);
        tree_cell_combined{end+1,1} = cell2vec(tree_cell);% get rid of all spaces and tabs 
        
%         max_width = max(size(tree_cell, 2), size(tree_header, 2));
%         if(max_width>size(tree_header, 2))
%             tree_header{end,max_width} = ''; 
%         end
%         if(max_width>size(tree_cell, 2))
%             tree_cell{end,max_width} = ''; 
%         end
%         tree_cell_combined = [tree_header(1:tree_line-1,:)' tree_cell']';         
        savecellfile(tree_cell_combined, cur_output_file); 

        species_cell_combined = repmat( get(AllSubTrees{i}, 'leafnames'), 1, 2); % New! save also file with species names !!! 
        for j=1:size(species_cell_combined,1)
            species_cell_combined{j,1} = genome_ver_to_organism(species_cell_combined{j,2}); 
        end
        savecellfile(species_cell_combined, strrep(cur_output_file, '.mod', '_species.txt')); % save species names 
        
    end            
end
    


