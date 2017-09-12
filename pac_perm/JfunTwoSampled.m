% This is just some vector function whose zeroes we want to find
% numerically. x_vec is here (X, t, s). f_star is treated as a parameter
% !!! 
% Try without the factor 2 !!!!
function [J , y_vec]  = JfunTwoSampled(x_vec, sigma, alpha, f_star, one_side_flag, soft_constrain_flag)

TOL = 0.0000000000001;

TWO = 2; % Should be two !!!!!. Just checking without !!!
I = sqrt(-1);

% Enough with the toy example, now compute the true 3-d function : 
C_alpha = norminv(1-0.5*alpha);   % Compute C_alpha from alpha

% Try new function with soft constrains !!! 

BIG_J=10;


if(soft_constrain_flag)
    minus_iy = log(x_vec(2)+1)/(2*BIG_J);
    minus_iz = (log(x_vec(3)+1)-log(x_vec(2)+1))/(2*BIG_J);
    minus_iw = minus_iz; % dummy - wrong !!! 
else
    minus_iy = 0;
    minus_iz = 0;
    minus_iw = 0; 
end 


% SaddleHelperGaussianQuadlManyVars(a, b, TOL, x1, x2, sigma, u, t, s, deriv_flag)
% y(1) = F_z/i y(2) = (F_y - F_z)/i y(3) = -F_x
% These are the important saddle equations : 
y_vec(1) = -TWO .* SaddleHelperGaussianQuadlManyVars(0, 99,TOL, x_vec(1), x_vec(2), sigma, x_vec(3), x_vec(4), x_vec(5), 1);    % F_x1
y_vec(2) = -TWO .* SaddleHelperGaussianQuadlManyVars(0, 99, TOL, x_vec(1), x_vec(2), sigma, x_vec(3), x_vec(4), x_vec(5), 2);   % F_x2
y_vec(3) = -TWO .* SaddleHelperGaussianQuadlManyVars(0, 99, TOL, x_vec(1), x_vec(2), sigma, x_vec(3),  x_vec(4), x_vec(5), 3) - (1-alpha)*I ; % F_y
y_vec(4) = -TWO .* SaddleHelperGaussianQuadlManyVars(0, 99, TOL, x_vec(1), x_vec(2), sigma, x_vec(3),  x_vec(4), x_vec(5), 4) - (1-alpha)*I; % F_z
y_vec(5) = -TWO .* SaddleHelperGaussianQuadlManyVars(0, 99, TOL, x_vec(1), x_vec(2), sigma, x_vec(3),  x_vec(4), x_vec(5), 5)- (1-alpha*f_star)*I; % F_w

% Here we try to see if the derivatives make sense
This_Should_Be_Zero_At_Saddle = y_vec; 

% Now divide according to two parts
%%%yv3 = [TWO .* SaddleHelperGaussianQuadlManyVars(0, 99, TOL, x_vec(1), x_vec(2), sigma, x_vec(3),  x_vec(4), x_vec(5), 3)  ]

y_vec = y_vec';

% Now compute the Jacobian. Note : Here we do differentiation with respect
% to t and s, just for numerical computation. This is different than doing 
% differentiation with respect to y and z and w !!! 
% This J is not symmetric !!! 
J = zeros(5);


for i=1:5
    for j=1:5
        J(i,j) = -2*SaddleHelperGaussianQuadlManyVars(0, 199, TOL, x_vec(1), x_vec(2), sigma, x_vec(3), x_vec(4), x_vec(5), 1000+10*i+j); % Here take integral from 0 to infinity
    end
end
