% Perform unique and count how many times each value appears
% 
% Input: 
% x - vector of values
% flag - additional flag for unique function (e.g. 'rows')
% 
% Output: 
% U - unique vector
% h - histogram of counts
% 
function [U h] = unique_with_counts(x, flag, varargin)

if(exist('flag', 'var'))
    [U I J] = unique(x, flag); 
else
    [U I J] = unique(x); 
end
h = hist(J, 1:length(U)); 

