% Get startionary distribution for a Makrov chain
% Input: 
% M - transition matrix of Markov chain
% Output: 
% pi - stationary distribution
% 
function pi = markov_chain_stationary_dist(M)
[pi, ~] = eigs(M', 1); pi = pi ./ sum(pi); 
