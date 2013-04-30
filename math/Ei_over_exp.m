% The standard exponential integral function divided by exp(x) (for asymptotic expansions)
% 
% Input: 
% x - variable
% 
% Output: 
% ret - equal to integral_{-infinity}^x  exp(t)/t dt / exp(x)
%
function ret = Ei_over_exp(x)

ret = Ei(x) / exp(x); 

if(isnan(ret))
    ret = 1/x;
end



