% Get startionary distribution for a Makrov chain
function pi = markov_chain_stationary_dist(M)
[pi, ~] = eigs(M', 1); pi = pi ./ sum(pi); 
