% Compute function given by a circuit (a general function - should move outside !!!)
% circuit has edges on non-diagonal elements and gate type on diagonal elements
% Problem: each gate also needs parameters (e.g. affine gate) - not just
% one type
%
% Input:
% x - input values for the first N variables
% circuit - a DAG with circuit representation. Nodes values represent
%           circuit type and edges represent inputs
% circuit_gate_params - cell array with parameters for each gate
% num_outputs - how many outputs to generate
%
% Output:
% z - output value
%
function  z = apply_circuit(x, circuit, circuit_gate_params, num_outputs)

[num_outputs N] = size(x); % number of input loci
m = length(circuit); % total number of vertices
V = zeros(num_outputs,m)-1; % initilize vector
V(:,1:N) = x; % get inputs

circuit_roots = get_graph_sources(circuit);
[node_dists node_discovery_time pred] = bfs(sparse(circuit),circuit_roots(1)); % determine order using breadth first search

[sorted_node_dists sort_perm] = sort(node_dists); % sort to get ordering (takes care of everything except unreachable nodes
first_reachable_node = find(sorted_node_dists > 0, 1); % take first reachable node (not the souce node)

for i=1:length(circuit_roots) % first evaluate root gates (no input gates)
    cur_node = circuit_roots(i); 
    if(circuit(cur_node,cur_node) > 0) % only roots which have a gate 
        cur_parents = setdiff(parents(circuit, cur_node), cur_node); % go to children
        V(:,cur_node) = apply_gate(V(:,cur_parents), circuit(cur_node,cur_node), ...
            circuit_gate_params{cur_node});   % apply gate
    end
end

for i=first_reachable_node:m % Traverse the circuit in bfs order - not the '-1' nodes (not visited yet)
    cur_node = sort_perm(i);
    if(min(V(:,cur_node)) < 0) % test if i'th node already assigned
        cur_parents = setdiff(parents(circuit, cur_node), cur_node); % go to children
        
        V(:,cur_node) = apply_gate(V(:,cur_parents), circuit(cur_node,cur_node), ...
            circuit_gate_params{cur_node});   % apply gate
    end
end

circuit_sink = get_graph_sinks(circuit);
z = V(:,circuit_sink); % output value at sink (supposed to be only one sink!)



