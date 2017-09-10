% Generate a uniform spanning tree on n vertices
% The appropriate thing to do is use Wilson's algorithm 
% with loop erased random walks, but here we give a simpler 
% implementation which gives you 'psuedo' uniform trees
% (so the actual sampling distribution might deviate from uniform)
function Tree = uniform_spanning_tree(n, directed)

Tree = zeros(n); Tree(1,2) = 1; % First edge
r = ceil(rand(1,n) .* [0:n-1]);
for i=3:n   
    Tree(r(i),i) = 1;     % Connect to a random previously connected edge    
end

P = randperm(n); % random ordering. This give more weight to self-isomorphic trees
Tree = Tree(P,P); 
if(directed == 0)
    Tree = Tree+Tree';
end






