% Compute variance for a histogram
% Input: 
% x - values
% p - frequencies
% 
% Output: 
% v - variance of distribution 
% 
function v = var_hist(x, p)

v =  max(0, mean_hist(x.^2, p) - mean_hist(x,p)^2);
