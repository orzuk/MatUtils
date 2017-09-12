% The standard exponential integral function
% 
% Input: 
% x - variable
% 
% Output: 
% ret - equal to integral_{-infinity}^x  exp(t)/t dt
%
function ret = Ei(x)

ret = (-expint(-x) - 1i*pi);
