% Returns the mean of a truncated chi square, E[max(W, c)] where W ~ chi2(k)
% Input: 
% c - value of truncation
% k - degrees of freedom
% 
function M = trunc_chi2stat(c, k)

M = c.*chi2cdf(c, k) + k.*(1-chi2cdf(c, k+2)); 
