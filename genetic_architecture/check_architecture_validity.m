% Test that the architecture gives always probabilities between zero and one
% 
% Input: 
% architecture_str - string representing architecture type
% params_struct - structure with architecture parameters
% N - number of pathways 
% 
% Output: 
% valid_flag - 1 if prob. is in [0,1] and 0 otherwise  
%
function valid_flag = check_architecture_validity(architecture_str, params_struct, N)

z_0 = genetic_architecture(zeros(1, N), architecture_str, params_struct, 1);
z_1 = genetic_architecture(ones(1, N), architecture_str, params_struct, 1);

valid_flag = ((z_1 <= 1) & (z_0 >= 0));



