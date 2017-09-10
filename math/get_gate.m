% Build a graph representing a computational gate
% 
% Input: 
% gate_type - type of gate (number)
% n - number of inputs 
% gate_params - ???
% 
% Output: 
% E - graph representing the gate
% 
function E = get_gate(gate_type, n, gate_params)

E = zeros(n+1); % last node is output
E(1:n,n+1) = 1; % get edges
E(n+1,n+1) = gate_type; 


