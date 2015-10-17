% Makes a cell array with identical copies of x
% 
% Input: 
% x - value
% m - width
% n - length
% 
% Output: 
% c - cell array of size mXn with x value replicated
%
function c = repcell(x, m, n)
c = repmat(x, m, n); 
c = mat2cell(c, ones(m,1)*size(c,1)/m, ones(n,1)*size(c,2)/n); 
