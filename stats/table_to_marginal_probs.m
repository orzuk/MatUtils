% Compute product of marginals for a probability table
% Used as expected values when assuming no associations between rows and columns
% Input: 
% T - matrix of joint probabilities
% Output: 
% T_marginal - matrix with marginals 
% 
function T_marginal = table_to_marginal_probs(T)

T_marginal = sum(T,2) * sum(T) ./ sum(T(:));  

