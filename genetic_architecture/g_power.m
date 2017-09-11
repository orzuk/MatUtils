% The function g(lambda) appearing in approximate power calculations 
function ret = g_power(lambda)

ret = (1+lambda).*log(1+lambda)-lambda;
