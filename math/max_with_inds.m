% Gets the maximum of two arrays, and also outputs indices 
% which were chose to be the maximal ones 
%
% The input: 
% a - first vec
% b - second vec
% 
% The output: 
% m - maximum vec
% inds - binary vec saying who was bigger
%
function [m inds] = max_with_inds(a, b)

m = max(a,b); 
inds = (m == b); % 1 means b is the max and zero means a is the max



