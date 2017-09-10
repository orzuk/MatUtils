% Get the maximum element in each member of a cell array
% 
% Input: 
% c - cell array
%  
% Output: 
% M - maximum value in each cell
% I - index of maximum value in each cell
% 
function [M, I] = max_cell(c)

n = length(c);
M = zeros(n,1); I=M;
for i=1:n
    [M(i), I(i)] = max(c{i}(:));
end
