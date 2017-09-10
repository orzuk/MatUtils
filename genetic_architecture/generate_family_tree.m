% Generate a family tree (full binary tree) with a given set of vertices
%
% Input:
% max_generation - depth of tree (how many generations)
%
% Output:
% family_tree - binary tree of family
%
function family_tree = generate_family_tree(max_generations) % New: compute sibling risk by simulating an entire family tree (should be different then risk for son!!)

AssignGeneralConstants;
N = 2^max_generations; % number grows exponentially
I = []; J = []; Vals = []; current_parents = 1:N;
for i=1:max_generations % go down level by level in the tree
    current_children = (current_parents(end)+1):(current_parents(end)+length(current_parents));  % / 2
    I = [I current_parents(1:2:end) current_parents(2:2:end) ...
        current_parents(1:2:end) current_parents(2:2:end)];
    J = [J current_children(1:2:end) current_children(1:2:end) ...
        current_children(2:2:end) current_children(2:2:end)];
    Vals = [Vals ones(1,2*length(current_children))];
    
    % Need to add also a twin at the last stage .. 
    current_parents = current_children(2:2:end);    
end
    unique_I = 1:max(I); % add gender information 
%    I = [I unique_I]; J = [J unique_I];
%    gender_vec = MALE + zeros(length(unique_I), 1); 
%    gender_vec(2:2:end) = FEMALE; Vals = [Vals gender_vec']; 

N = max(max(I), max(J));

family_tree = sparse(I, J, Vals, N, N); 


