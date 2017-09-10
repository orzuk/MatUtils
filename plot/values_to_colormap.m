% Compute colormap values for a set of indices 
% Input: 
% x - vector (matrix?) of real values 
% 
% Ouput: 
% c - Nx3 matrix of colors 
% w - vector of line weights / node sizes 
%
function [c w] = values_to_colormap(x, colormap_str, flip_flag)

if(~exist('colormap_str', 'var') || isempty(colormap_str))
    colormap_str = 'jet';
end
ag = findobj; num_open_figures = max(ag(find(ag==fix(ag)))); % Get number of open figures 
colormap(colormap_str);
 %
if(flip_flag)
    colormap(flipud(colormap));
end
% colormap('default'); 
color_mat = colormap; ncolors = size(color_mat,1); % get 64 colors
if(num_open_figures==0)
    close all;
end
    
    
min_x = min(x(:)); max_x = max(x(:)); range_x = max_x - min_x; % Compute range. why take only positive values? min(x(x>0))


c = color_mat(max(1,ceil( ncolors * (x - min_x) / range_x  )), :); % 
w = x ./ min_x;



