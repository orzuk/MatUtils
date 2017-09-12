% Set number of figures in each row and column in subplot 
% Input: 
% n - total number of figures
% 
% Output: 
% num_h - number in each row
% num_w - number in each column
% 
function [num_h, num_w] = num_to_height_width(n)

num_h = ceil(sqrt(n)); num_w = floor(sqrt(n));
if(num_h*num_w < n)
    num_w = num_h;
end