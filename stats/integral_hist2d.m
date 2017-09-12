% Perform two-dimensional integral for a histogram (we use simple rectangle integration)
%
% Input:
% x - values on x-axis
% y - values on y-axis
% p - their density (matrix)
%
% Output:
% s - the integral  \sum_{i,j} p(i,j) * \delta_x * \delta_y
%
function s = integral_hist2d(x, y, p)

% s = 0.5 .* sum(diff(x) .* (p(1:end-1) +  p(2:end)));

if(~isvector(x)) % Here x and y are already given as a mesh-grid
    x = x(1,:); y = y(:,1); % take one representative vector
end
s = 0.25 .* sum( sum( (vec2column(diff(x)) * vec2row(diff(y))) .* ...
    (p(1:end-1,1:end-1) + p(2:end,1:end-1) + p(1:end-1,2:end) + p(2:end,2:end)) ) );
