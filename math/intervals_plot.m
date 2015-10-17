% Plot a set of intervals horizontaly on different heights
% 
% Input: 
% start_pos - interval starts
% end_pos - interval ends 
% hegiht - height of each interval 
% color_str - interval colors
% 
function Dummy = intervals_plot(start_pos, end_pos, height, color_str, line_width, varargin)

if(~exist('color_str', 'var')) % set default color 
    color_str = 'b';
end
if(~exist('line_width', 'var') || isempty(line_width))
    line_width = 1;
end
n = length(start_pos); hold on;

% first plot the end-points
plot(vec2column(start_pos), vec2column(height) + zeros(n,1), ['*', color_str]);
plot(vec2column(end_pos),  vec2column(height) + zeros(n,1), ['*', color_str]);

% Now plot the intervals themselves
if(length(height) == 1) % one height
    height = repmat(height, 2, n);
else % n different heights
    height = repmat(vec2row(height), 2, 1); 
end
Dummy = 0; line([vec2column(start_pos) vec2column(end_pos) ]', height + zeros(2, n), ...
    'color', color_str, 'linewidth', line_width); 

