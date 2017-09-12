% Perform an integral by deviding into small areas to increase accuracy
function F = SaddleHelperGaussianQuadl(a, b, TOL, x, sigma, t, deriv_flag)

% Now try to seperate the integral into 3
c = a + (b-a)/100; % Note that a shouldn't be negative !!! 
d = max(a+(b-a)/10, c);
F = quadl('SaddleHelperGaussian', a, c, TOL, [], x, sigma, t, deriv_flag) + ...
    quadl('SaddleHelperGaussian', c, d, TOL, [], x, sigma, t, deriv_flag) + ...
    quadl('SaddleHelperGaussian', d, b, TOL, [], x, sigma, t, deriv_flag);


