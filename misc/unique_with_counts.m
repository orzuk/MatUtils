% Perform unique and count how many times each value appears
% 
% Input: 
% x - vector of values
% flag - additional flag for unique function (e.g. 'rows')
% weights - NEW! allow weigthed counts. Vector of same size as x with weights 
% 
% Output: 
% U - unique vector
% h - histogram of counts
% 
function [U, h] = unique_with_counts(x, flag, weights, varargin)

if(exist('flag', 'var') && (~isempty(flag)))
    [U, ~, J] = unique(x, flag); 
else
    [U, ~, J] = unique(x); 
end
if(~exist('weights', 'var') || isempty(weights))
    h = hist(J, 1:length(U)); 
else
    h = weighted_hist(J, weights, 1:length(U)); 
end
if(iscolumn(U)) % set same dimension
    h = vec2column(h);
end
