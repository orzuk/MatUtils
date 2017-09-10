% Compute the running/cummulative min of a vector, 
% by simply calling the cummulative max of the matlab utilities package 
% Input: 
% x - vector
% dim - dimension (optional) 
% 
% Output: 
% y - vector of cumulative minimum values 
function y = cummin(x, dim, varargin)

if (nargin < 2)    
    y = -cummax(-x); % call the cummax function 
else
    y = -cummax(-x, dim); 
end
    
    
    
    
   