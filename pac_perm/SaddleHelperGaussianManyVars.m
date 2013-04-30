% This shows the weird function that we should integrate in the saddle point equations. 
% The derive flag controls between the many kinds of densities.
function F = SaddleHelperGaussianManyVars( c, x_vec, sigma, deriv_flag)

%%%%%Inside_SaddleHelperGaussianManyVars_sigma = sigma

% In the function we need x_vec so that we have few variables in order not
% to upset the quadl function.
x1 = x_vec(1); x2 = x_vec(2); u = x_vec(3); t = x_vec(4); s = x_vec(5);

I = sqrt(-1);

q_c =  (1./sqrt(2.*pi)) .* exp(-c.^2./2); % Set the density q to be gaussian

% Function for F :
if(deriv_flag == 0)
    F = q_c .* log(A(x1, x2, sigma, u, t, s, c));
end


% Five functions for saddle point : 
if(deriv_flag == 1)
    F = q_c .* Px(x1,c, sigma) .* ( PInt(x2,c, sigma) .* (u-t-s)+t) ./ A(x1,x2,sigma,u,t,s,c);
end
if(deriv_flag == 2)
    F = q_c .* Px(x2,c, sigma) .* ( PInt(x1,c, sigma) .* (u-t-s)+s) ./ A(x1,x2,sigma,u,t,s,c);
end
if(deriv_flag == 3)
    F = (-I) .* q_c .* PInt(x1,c, sigma) .* ( PInt(x2,c, sigma) .* (u-t)+t+1) ./ A(x1,x2,sigma,u,t,s,c);
end
if(deriv_flag == 4)
    F = (-I) .* q_c .* PInt(x2,c, sigma) .* ( PInt(x1,c, sigma) .* (u-s)+s+1) ./ A(x1,x2,sigma,u,t,s,c);
end
if(deriv_flag == 5)
    F = (-I) .* q_c .* [PInt(x1,c, sigma) .* PInt(x2,c, sigma) .* (u-t-s-1) + PInt(x1,c, sigma) .* (t+1) + PInt(x2,c, sigma) .* (s+1)] ./ A(x1,x2,sigma,u,t,s,c);
end


% Special functions for two samples - here we must have s=t=u=0, and
% x_1=x_2=x_alpha !!!
if((deriv_flag > 80) && (deriv_flag < 90))% Here deal with the gradient of F
    var_flag = mod(deriv_flag, 10);  % Choose variable
    F = q_c .* A_1st_deriv(x1,x2,sigma,u,t,s,var_flag,c) ./ A(x1,x2,sigma,u,t,s, c);        % (x_1,x_2,sigma,u,t,s, var_flag, c)
end

if((deriv_flag > 899) && (deriv_flag < 1000))% Here deal with the gradient of F
    var_flag1 = mod((deriv_flag-mod(deriv_flag,10))/10,10);  % Choose 1st variable
    var_flag2 = mod(deriv_flag, 10);  % Choose 2nd variable
    F = q_c .* (A_2nd_deriv(x1, x2, sigma, u, t, s, var_flag1, var_flag2, c) .*  A(x1,x2,sigma,u,t,s, c) - ...
        A_1st_deriv(x1,x2,sigma,u,t,s,var_flag1,c) .* A_1st_deriv(x1,x2,sigma,u,t,s,var_flag2,c)) ./ (A(x1,x2,sigma,u,t,s, c) .^ 2); % Not dividing !!! 
end

% Here do derivatives with respect to u,t,s for the Newton method
if((deriv_flag > 999) && (deriv_flag < 1100))% Here deal with the gradient of F
    var_flag1 = mod((deriv_flag-mod(deriv_flag,10))/10,10);  % Choose 1st variable
    var_flag2 = mod(deriv_flag, 10);  % Choose 2nd variable
    F = q_c .* (A_2nd_deriv_for_Newton(x1, x2, sigma, u, t, s, var_flag1, var_flag2, c) .*  A(x1,x2,sigma,u,t,s, c) - ...
        A_1st_deriv(x1,x2,sigma,u,t,s,var_flag1,c) .* A_1st_deriv_for_Newton(x1,x2,sigma,u,t,s,var_flag2,c)) ./ (A(x1,x2,sigma,u,t,s, c) .^ 2);
    % In the above line, the first A_deriv is to y,z,w and the second
    % is to u,t,s

end


