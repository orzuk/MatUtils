% Compute skewness for a histogram
%
% Input: 
% x - values
% p - their probabilities
% 
% Output: 
% mu - mean value 
% 
function s = skewness_hist(x, p)

s = moment_hist(x, p, 3) / std_hist(x, p).^3; 
