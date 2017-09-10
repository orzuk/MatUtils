% Generate random points inside a 3-d sphere
% Input: 
% n - number of points
% R - radius (optional) 
%
% Output: 
% x - vector of points 
% 
function x = randsphere(n, R, varargin)

x = randn(3,n);
d = sqrt( x(1,:).^2+x(2,:).^2+x(3,:).^2 );
x(1,:) = x(1,:)./d;
x(2,:) = x(2,:)./d;
x(3,:) = x(3,:)./d;

r = (rand(1,n)).^(1/3);
x = x .* repmat(r, 3, 1);

if(exist('R', 'var'))
    x = x .* R;
end
