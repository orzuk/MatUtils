% Build a directed graph represening the architecture
% 
% Input: 
% N - number of input variables
% architecture_str - type of architecture
% params_struct - architecture parameters
% 
% Output: 
% E - represents the edges of the (directed) circuit
% circuit_gate_params - parameters of each gate in the circuit
% 
function [E circuit_gates_params] = ...
    set_architecture_circuit(N, architecture_str, params_struct)

num_extra_nodes = sum(params_struct.num_clasues); % intermidiate and output nodes

E = zeros(N+num_extra_nodes);
circuit_gates_params = cell(N+num_extra_nodes,1); % parameters for each gate

ctr = N+1;
cur_level_inds = 1:N;
for i=1:params_struct.levels % how 'deep' is the architecture
    new_level_inds = [];
    for j=1:params_struct.num_clasues(i)
        cur_inds = ((j-1)*params_struct.k_in_clause(i)+1: j*params_struct.k_in_clause(i));
        cur_inds = cur_level_inds(cur_inds);
        input_gates = diag(E(cur_inds,cur_inds));
        E([cur_inds ctr], [cur_inds ctr]) = E([cur_inds ctr], [cur_inds ctr]) + ...
            get_gate(params_struct.gate_type(i), ...
            params_struct.k_in_clause(i), params_struct.gate_params(i));
        %diag(E(cur_inds,cur_inds)) = input_gates; % return gates of input vertices to current gate
        circuit_gates_params{ctr}.linear_coef_vec = params_struct.linear_coef_vec;
        if(isfield(params_struct, 'a'))
            circuit_gates_params{ctr}.a = params_struct.a(i);
        end
        if(isfield(params_struct, 'b'))
            circuit_gates_params{ctr}.b = params_struct.b(i);
        end
        new_level_inds = [new_level_inds ctr]; % build indices for next level
        ctr=ctr+1;
    end
    cur_level_inds = new_level_inds; % update inputs for next level
end

