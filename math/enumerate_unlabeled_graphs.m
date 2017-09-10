% Enumerate all non-isomorphic undirected graphs 
function non_isomorphic_graphs = enumerate_unlabeled_graphs(N)

all_perms = perms(1:N); % Take all permutations

num_edges = N*(N-1)/2;

graph_perms = zeros(factorial(N), num_edges); % permutation on edges

ind=1;
% loop on edges
for i=1:factorial(N)
    cur_edge=1;
    for j=1:N
        for k=j+1:N
            permed_j = min(all_perms(i,j), all_perms(i,k));
            permed_k = max(all_perms(i,j), all_perms(i,k));

            % get the lexicographic place of the (i,j) edge
            graph_perms(i,cur_edge) = (permed_j-1)*(2*N-permed_j)/2 + permed_k-permed_j; cur_edge=cur_edge+1;
        end
    end
end

% Now generate all graphs
all_graphs = zeros(2^num_edges, num_edges);

for i=1:num_edges
    all_graphs(:,i) = mod( floor([0:2^num_edges-1] ./ (2^(i-1))), 2);
end


% Now start assigning colors
colors = zeros(2^num_edges,1);

% Go over all graphs
color=1;
for i=1:2^num_edges
    if(colors(i) == 0) % go only on uncolored graphs
        cur_graph = all_graphs(i,:);
        all_graph_perms = cur_graph(graph_perms);
       
        all_graph_indices = sum(all_graph_perms .* repmat(2.^[0:num_edges-1],factorial(N), 1), 2) + 1;

        colors(all_graph_indices)=color; color=color+1;
    end
end

[B,I,J]=unique(colors);
non_isomorphic_graphs = all_graphs(I,:);
size(non_isomorphic_graphs)
