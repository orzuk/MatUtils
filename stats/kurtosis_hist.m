% Compute kurtosis for a histogram
%
% Input: 
% x - values
% p - their probabilities
% 
% Output: 
% mu - mean value 
% 
function k = kurtosis_hist(x, p)

k = moment_hist(x, p, 4) / std_hist(x, p).^4 - 3; 
