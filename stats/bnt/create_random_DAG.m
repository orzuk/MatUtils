function dag = create_random_DAG(nodes, edge_prob)
% create_random_DAG creates a 'random' directed acyclic graph.
% The number of input nodes is nodes.
% For each pair of nodes (i,j), there will be an edge with probability
% edge_prob. This is : edge_prob/2 for i->j and edge_prob/2 for j->i
% Method : 
% 1. Randomize a permutation
% 2. Flip a coin for every edge, in the direction of this permutation
%

% create order permutation
p = randperm(nodes);

dag = zeros(nodes);
% flip coins for edges
tmp = rand(nodes) < edge_prob;

dag(find(tmp)) = 1;


% remove lower triangle
x = [1:nodes*nodes];
dag(find(floor( (x-1)/nodes) <= rem(x-1, nodes))) = 0;

% arrange according to the perm
dag = dag(p, p);
