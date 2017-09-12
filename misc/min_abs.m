% Trancate x at M in absolute value
function x_min = min_abs(x, m)

x_min =  max(min(x, m), -m); 

