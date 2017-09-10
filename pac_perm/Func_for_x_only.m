% Compute saddle point when changing only x
function F = FUNC_for_x_only( x, f_star, alpha, C_alpha, sigma)
TOL = 0.0000000001;


F = 2 .* quadl('SaddleHelperGaussian', C_alpha, 999, TOL, [], x, sigma, 0, 0) - alpha .* (1-f_star);