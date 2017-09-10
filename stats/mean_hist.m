% Compute mean for a histogram
% 
% Input: 
% x - values
% p - probabilities
% 
% Output: 
% mu - mean of distribution 
% 
function mu = mean_hist(x, p)

mu = sum(x.*p) / sum(p); 
