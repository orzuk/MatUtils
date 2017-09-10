% Find value with highest count in a vector.
% This is like matlab's mode function 
% 
% Input: 
% x - vector of different values
% 
% Output: 
% m - most common value
% count - how many times it appears
% 
function [m count]= majority(x)
[u h] = unique_with_counts(x);
[count I] = max(h);
m = u(I);

