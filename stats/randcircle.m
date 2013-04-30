% Draw random points in a circle
% 
% Input:
% n - number of points 
% R - radius (optional)
%
% Output: 
% x - vector of points 
%
function x = randcircle(n, R, varargin)

% x = randn(2,n);
% d = sqrt( x(1,:).^2+x(2,:).^2 );
% x(1,:) = x(1,:)./d;
% x(2,:) = x(2,:)./d;
% 
% if(exist('R', 'var'))
%     x = x .* R;
% end

r = sqrt(rand(n,1)); 
theta = 2*pi.*rand(n, 1); % randomize angle

x = zeros(2,n);
x(1,:) = r .* cos(theta);
x(2,:) = r .* sin(theta);



