% Take modulus but if it's zero set it as n (instead of zero)
% 
% Input: 
% x - value
% n - modulus
% 
% Output: 
% y - (x-1) mod n + 1
%
function y = mod_max(x, n)
y = mod(x,n); 
if(y == 0) 
    y = n;
end

