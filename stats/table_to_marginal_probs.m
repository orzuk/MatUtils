% Compute product of marginals for a probability table
% Used as expected values when assuming no associations between rows and columns
function T_marginal = table_to_marginal_probs(T)

T_marginal = sum(T,2) * sum(T) ./ sum(T(:));  

