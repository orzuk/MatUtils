% Generate an adjancancy matrix for a full balanced binary tree
function E = binary_tree(depth)

num_nodes = 2^depth-1; num_edges = 2^depth-2; 
ones_vec = ones(num_edges,1); 

I = mat2vec(repmat(1:(num_edges/2), 2, 1)); 
J = I.*2; J(2:2:end) = J(2:2:end)+1; 

E = sparse(I, J, ones_vec, num_nodes, num_nodes);



