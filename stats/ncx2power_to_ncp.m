% Compute the non-centrality parameter needed to achieve a desired power
% for non-central chi-square distribution
% 
% Input: 
% alpha - significance level
% v - degrees of freedom
% pow - desired power
% 
% Output
% NCP - non-centrality parameter achieving desired power 
% 
function NCP = ncx2power_to_ncp(alpha, v, pow)

x_alpha = chi2inv(1-alpha, v); % compute threshold 

options = optimset('tolx', 0.0000000000000000000001); 
NCP = fminbnd(@(x_ncp) abs(ncx2cdf(x_alpha, v, x_ncp)+pow-1), 0, 3*x_alpha+5, options);

