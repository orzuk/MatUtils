% This shows the joint probability distribution. It depends on the original
% correlation function. 
function F = SaddleHelperIntGaussian(x, sigma, t,  C_alpha, one_side_flag)


F = quad('SaddleHelperGaussian', C_alpha, 999, [], [],  x,  sigma, t);
