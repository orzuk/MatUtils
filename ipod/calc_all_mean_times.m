% Calculate mean time until a fraction alpha of the popultaion is infected
% starting from ANY initial state. Computation currently
% is exponential, and thus we must restrict ourselves to small values of N
%
% Input: 
% E - graph adjacency matrix  (sparse)
% X_init - initial state of all nodes
% alpha - desired fraction of infected nodes
% lambda - ?
%
% Output: 
% t_mean - mean time (# steps) until at least alpha of nodes is infected
% 
function t_mean =calc_all_mean_times(E, X_init, alpha, lambda)

N = length(E) % number of vertices

t_mean = zeros(1,2^N); % the mean time
ham_vec = hamming_weight([0:2^N-1]);  % Hold the hamming vec of each subset
phi_edges_vec = zeros(1,2^N); % For each state store the edges with one vertex in the set and one in its complement
phi_vertex_vec = zeros(1,2^N); % For each state store the vertices which are adjacent to it (this is in fact the boundary of the complement)


% Start the recursion
for num_on=ceil(alpha*N)-1:-1:1  % loop over the number of on vertices in backwards order
    cur_states_indexes = find(ham_vec == num_on); % Find the relevant states for this stage

    t_mean(cur_states_indexes) =  1/lambda;

    for i=1:N     % loop over all possible neighbors, assume for now that N <=32
        add_i_states_indexes = bitor(cur_states_indexes-1, 2^(i-1))+1; % The indexes of the sets with one more vertex

        % Compute the degree-into the set A
        i_neighbors = sum(2.^(find(E(i,:))-1));
        d_i_indexes = hamming_weight(bitand(i_neighbors, cur_states_indexes-1)) .* (1- mod( floor((cur_states_indexes-1) ./ 2^(i-1)), 2));

        t_mean(cur_states_indexes) = t_mean(cur_states_indexes)  + d_i_indexes .* t_mean(add_i_states_indexes);
        phi_edges_vec(cur_states_indexes) = phi_edges_vec(cur_states_indexes) + d_i_indexes; % add the number of edges to the count
    end

    num_on_is = num_on

    % Divide by the boundary size
    t_mean(cur_states_indexes) = t_mean(cur_states_indexes) ./ phi_edges_vec(cur_states_indexes);

end  % loop on steps








