% Perform an integral by deviding into small areas to increase accuracy
function F = SaddleHelperGaussianQuadlManyVars(a, b, TOL, x1, x2, sigma, u, t, s, deriv_flag)
             


x_vec = zeros(1,5);
x_vec(1) = x1; x_vec(2) = x2; x_vec(3) = u; x_vec(4) = t; x_vec(5) = s;







% Now try to seperate the integral into 3
% c = a + (b-a)/100; % Note that a shouldn't be negative !!! 
% d = max(a+(b-a)/10, c);
% F = quadl('SaddleHelperGaussianManyVars', a, c, TOL, [], x_vec, sigma, deriv_flag) + ...
%     quadl('SaddleHelperGaussianManyVars', c, d, TOL, [], x_vec, sigma, deriv_flag) + ...
%     quadl('SaddleHelperGaussianManyVars', d, b, TOL, [], x_vec, sigma, deriv_flag);




% Try to use many sub-intervals
subs_num = 10;
dum_a = 1; dum_b = 1 + b-a;
subs = 0:(log(dum_b)/(subs_num-1)):log(dum_b);
subs = exp(subs)+a-1;

F = 0;
for i=1:subs_num-1
    F = F + quadl('SaddleHelperGaussianManyVars', subs(i), subs(i+1), TOL, [], x_vec, sigma, deriv_flag);
end
