% This is just some vector function whose zeroes we want to find
% numerically. x_vec is here (X, t, s). f_star is treated as a parameter
% !!! 
% Try without the factor 2 !!!!
function [J , y_vec]  = Jfun(x_vec, sigma, alpha, f_star, one_side_flag, soft_constrain_flag)

% sigma
% alpha 
% f_star
% one_side_flag
% soft_constrain_flag

% Compute the value of the function 
%%% y_vec =  [sin(x_vec(1)+x_vec(3)), cos(x_vec(1)*x_vec(3))]';

% Compute the Jacobian
%%% J = [ cos(x_vec(1)+x_vec(3)) -x_vec(3) * sin(x_vec(1)*x_vec(3)); cos(x_vec(1)+x_vec(3)) -x_vec(1) * sin(x_vec(1)*x_vec(3)) ]; 
TOL = 0.0000000000001;

TWO = 2; % Should be two !!!!!. Just checking without !!!


% Enough with the toy example, now compute the true 3-d function : 
% alpha .* (1-f_star)
C_alpha = norminv(1-0.5*alpha);   % Compute C_alpha from alpha
% x is x_vec(1), t is x_vec(3), s is x_vec(2)

% Try new function with soft constrains !!! 

BIG_J=10;

% % % y_vec(1) = 2.* quadl('SaddleHelperGaussian', C_alpha, 99,TOL, [], x_vec(1), sigma, x_vec(2), 0) - alpha .* (1-f_star); 
% % % y_vec(2) = 2.* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 0) + alpha .* (1-f_star) - (1-alpha);  % here we take s !  
% % % y_vec(3) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 1) + ...
% % %            quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 1);  % here we took Px instead of PInt
% % % y_vec = y_vec';
% % % 
% % % % Now compute the Jacobian. Note : Here we do differentiation with respect
% % % % to t and s, just for numerical computation. This is different than doing 
% % % % differentiation with respect to y and z !!! 
% % % J = zeros(3);
% % % 
% % % J(1,1) = 2 .* quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 2);
% % % J(1,2) = 2 .* quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 3);
% % % J(1,3) = 0;
% % % 
% % % J(2,1) = 2 .* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 2);
% % % J(2,2) = 0;
% % % J(2,3) = 2 .* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 3);
% % % 
% % % J(3,1) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 4) + ...
% % %          quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 4);
% % % J(3,2) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 5);
% % % J(3,3) = quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 5);  % Here no factor of two should appear !!! 
% % % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





if(soft_constrain_flag)
    minus_iy = log(x_vec(2)+1)/(2*BIG_J);
    minus_iz = (log(x_vec(3)+1)-log(x_vec(2)+1))/(2*BIG_J);
else
    minus_iy = 0;
    minus_iz = 0;
end 

% % % % % % y_vec(1) = TWO.* quadl('SaddleHelperGaussian', C_alpha, 99,TOL, [], x_vec(1), sigma, x_vec(2), 0) - alpha .* (1-f_star); 
% % % % % % y_vec(2) = TWO.* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 0) + alpha .* (1-f_star) - ...
% % % % % %             (1-alpha);  % here we take s !  
% % % % % % y_vec(3) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 1) + ...
% % % % % %            quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 1);  % here we took Px instead of PInt
% % % % % % y_vec = y_vec';
% % % % % % 
% % % % % % % Now compute the Jacobian. Note : Here we do differentiation with respect
% % % % % % % to t and s, just for numerical computation. This is different than doing 
% % % % % % % differentiation with respect to y and z !!! 
% % % % % % J = zeros(3);
% % % % % % 
% % % % % % J(1,1) = TWO .* quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 2);
% % % % % % J(1,2) = TWO .* quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 3);
% % % % % % J(1,3) = 0;
% % % % % % 
% % % % % % J(2,1) = TWO .* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 2);
% % % % % % J(2,2) = 0;
% % % % % % J(2,3) = TWO .* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 3);
% % % % % % 
% % % % % % J(3,1) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 4) + ...
% % % % % %          quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 4);
% % % % % % J(3,2) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 5);
% % % % % % J(3,3) = quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 5);  % Here no factor of two should appear !!! 


% y(1) = F_z/i y(2) = (F_y - F_z)/i y(3) = -F_x
% These are the important saddle equations : 
y_vec(1) = TWO.* SaddleHelperGaussianQuadl(C_alpha, 99,TOL, x_vec(1), sigma, x_vec(3), 0) - alpha .* (1-f_star) + soft_constrain_flag .* minus_iz; 
y_vec(2) = TWO.* SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, x_vec(2), 0) + alpha .* (1-f_star) - (1-alpha) + ...
    soft_constrain_flag .* minus_iy - soft_constrain_flag .* minus_iz;  % here we take s !  
y_vec(3) = SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, x_vec(3), 1) + ...
           SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, x_vec(2), 1);  % here we took Px instead of PInt
y_vec = y_vec';

% Now compute the Jacobian. Note : Here we do differentiation with respect
% to t and s, just for numerical computation. This is different than doing 
% differentiation with respect to y and z !!! 
J = zeros(3);

J(1,1) = TWO .* SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, x_vec(3), 2);
J(1,2) = soft_constrain_flag .* 1./(2*(x_vec(3)+1)*BIG_J);
J(1,3) = TWO .* SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, x_vec(3), 3) -soft_constrain_flag .* 1/(2*(x_vec(2)+1)*BIG_J);

J(2,1) = TWO .* SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, x_vec(2), 2);
J(2,2) = TWO .* SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, x_vec(2), 3) -soft_constrain_flag .* 1/(2*(x_vec(3)+1)*BIG_J);
J(2,3) = 0;

J(3,1) = SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, x_vec(3), 4) + ...
         SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, x_vec(2), 4);
J(3,2) = SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, x_vec(2), 5);  % Here no factor of two should appear !!! 
J(3,3) = SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, x_vec(3), 5);



%J_IS = J;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% NEW IDea : Compute the Jacobian Stupidly using the Derivative's
% Definition !!!!!
% Stupid_J = zeros(3);
% resres = 0.0000001;
% 
% Stupid_y_vec(1) = TWO.* quadl('SaddleHelperGaussian', C_alpha, 99,TOL, [], x_vec(1)+resres, sigma, x_vec(2), 0) - alpha .* (1-f_star); 
% Stupid_y_vec(2) = TWO.* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1)+resres, sigma, x_vec(3), 0) + alpha .* (1-f_star) - (1-alpha);  % here we take s !  
% Stupid_y_vec(3) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1)+resres, sigma, x_vec(2), 1) + ...
%            quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1)+resres, sigma, x_vec(3), 1);  % here we took Px instead of PInt
% Stupid_J(1,1) = Stupid_y_vec(1)-y_vec(1);
% Stupid_J(1,2) = Stupid_y_vec(2)-y_vec(2);
% Stupid_J(1,3) = Stupid_y_vec(3)-y_vec(3);
%        
% 
% Stupid_y_vec(1) = TWO.* quadl('SaddleHelperGaussian', C_alpha, 99,TOL, [], x_vec(1), sigma, x_vec(2)+resres, 0) - alpha .* (1-f_star); 
% Stupid_y_vec(2) = TWO.* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 0) + alpha .* (1-f_star) - (1-alpha);  % here we take s !  
% Stupid_y_vec(3) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2)+resres, 1) + ...
%            quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3), 1);  % here we took Px instead of PInt
% Stupid_J(2,1) = Stupid_y_vec(1)-y_vec(1);
% Stupid_J(2,2) = Stupid_y_vec(2)-y_vec(2);
% Stupid_J(2,3) = Stupid_y_vec(3)-y_vec(3);
% 
% 
% Stupid_y_vec(1) = TWO.* quadl('SaddleHelperGaussian', C_alpha, 99,TOL, [], x_vec(1), sigma, x_vec(2), 0) - alpha .* (1-f_star); 
% Stupid_y_vec(2) = TWO.* quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3)+resres, 0) + alpha .* (1-f_star) - (1-alpha);  % here we take s !  
% Stupid_y_vec(3) = quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 1) + ...
%            quadl('SaddleHelperGaussian', 0, C_alpha, TOL, [], x_vec(1), sigma, x_vec(3)+resres, 1);  % here we took Px instead of PInt
% Stupid_J(3,1) = Stupid_y_vec(1)-y_vec(1);
% Stupid_J(3,2) = Stupid_y_vec(2)-y_vec(2);
% Stupid_J(3,3) = Stupid_y_vec(3)-y_vec(3);
% 
% Stupid_J = Stupid_J' ./ resres
% 
% J
           
%%%J=J'; % try .... ???? 