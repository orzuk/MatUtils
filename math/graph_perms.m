% Enumerate all graphs permutations on N nodes
function GP = graph_perms(N);

all_perms = perms(1:N); % Take all permutations

num_edges = N*(N-1)/2;

GP = zeros(factorial(N), num_edges); % permutation on edges

ind=1;
% loop on edges
for i=1:factorial(N)
    cur_edge=1;
    for j=1:N
        for k=j+1:N

            permed_j = min(all_perms(i,j), all_perms(i,k));
            permed_k = max(all_perms(i,j), all_perms(i,k));

            % get the lexicographic place of the (i,j) edge
            GP(i,cur_edge) = (permed_j-1)*(2*N-permed_j)/2 + permed_k-permed_j; cur_edge=cur_edge+1;

        end
    end
end



