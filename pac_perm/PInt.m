% This shows the joint probability distribution. It depends on the original
% correlation function. 
function F = PInt(x, c, sigma)

F = normcdf((x-c)/sigma) - normcdf((-x-c)/sigma);

