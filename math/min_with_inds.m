% Get the minimum of two arrays, and also output
% the indexes which were chose to be the minimal ones 
%
% Input: 
% a - first vec
% b - second vec
% 
% Output: 
% m - minimum vec
% inds - binary vec saying who was bigger
%
function [m inds] = min_with_inds(a, b)

m = min(a,b); 
inds = (m == b); % 1 means b is the min and zero means a is the min



