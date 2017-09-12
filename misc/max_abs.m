% Trancate x at M in absolute value
function x_max = max_abs(x, m)

x_max =  sign(x) .* max(abs(x), m);


