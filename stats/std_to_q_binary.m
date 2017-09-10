% Compute Bernoulli probability for a given standard deviation 
% 
% Input: 
% sigma - standard deviation
% 
% Output:
% q - binary one prob.
% 
function q = std_to_q_binary(sigma)

q  = (1-sqrt(1-4*sigma.^2))./2; % std of a binary variable
