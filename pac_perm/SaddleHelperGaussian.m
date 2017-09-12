% This shows the weird function that we should integrate in the saddle point equations. 
% The derive flag controls between the many kinds of densities.
function F = SaddleHelperGaussian( c, x,  sigma, t, deriv_flag)

q_c =  (1./sqrt(2.*pi)) .* exp(-c.^2./2); % Set the density q to be gaussian

% only two-sided flag so far. Our density q is standard gaussian
if(deriv_flag == 0)
    F = ( PInt(x,c,sigma) .* (t+1) .* q_c ) ./ (1 + PInt(x,c,sigma) .* t);
end
if(deriv_flag == 1)   % Here take Px
    F = ( Px(x,c,sigma) .* t .* q_c ) ./ (1 + PInt(x,c,sigma) .* t);
end

% Here we start with 2nd derivative stuff - The Jacobian : 
if(deriv_flag == 2) % here take Px
    F = ( Px(x,c,sigma) .* (t+1) .* q_c ) ./ (1 + PInt(x,c,sigma) .* t).^2;
end
if(deriv_flag == 3) % here take P*(1-P)
    F = ( PInt(x,c,sigma) .* (1 - PInt(x,c,sigma)) .* q_c ) ./ (1 + PInt(x,c,sigma) .* t).^2;
end
if(deriv_flag == 4) % here take complicated !!! 
    F = ( (Pxx(x,c,sigma) .* (1 + PInt(x,c,sigma) .* t) - t .* Px(x,c,sigma).^2).* t .* q_c ) ...
        ./ (1 + PInt(x,c,sigma) .* t).^2;
end
if(deriv_flag == 5) % here take Px without t+1
    F = ( Px(x,c,sigma) .* q_c ) ./ (1 + PInt(x,c,sigma) .* t).^2;
end

% This is a special flag for computing g. Log and no denominator ! 
if(deriv_flag == 6) 
    F = ( log(1 + PInt(x,c,sigma) .* t) .* q_c );
end

% Another special flag for the Jacobian of F Px (1-P)
if(deriv_flag == 7) % here take Px
    F = ( Px(x,c,sigma) .* (t+1) .* (1 - PInt(x,c,sigma)) .* q_c ) ./ (1 + PInt(x,c,sigma) .* t).^2;
end

% Special functions for two samples - here we must have s=t=u=0, and
% x_1=x_2=x_alpha !!!
if((deriv_flag > 80) && (deriv_flag < 90))% Here deal with the gradient of F
    var_flag = mod(deriv_flag, 10);  % Choose variable
    F = q_c .* A_1st_deriv(x,x,sigma,0,0,0,var_flag,c); %%% ./ A(x,x,sigma,0,0,0, c);        % (x_1,x_2,sigma,u,t,s, var_flag, c)
end

if((deriv_flag > 899) && (deriv_flag < 1000))% Here deal with the gradient of F
    var_flag1 = mod((deriv_flag-mod(deriv_flag,10))/10,10);  % Choose 1st variable
    var_flag2 = mod(deriv_flag, 10);  % Choose 2nd variable
%    F = (q_c .* (A_2nd_deriv(x, x, sigma, 0, 0, 0, var_flag1, var_flag2, c) .* A(x,x,sigma,0,0,0, c) - ...          %(x_1,x_2,sigma,u,t,s, var_flag1, var_flag2, c)
%        A_1st_deriv(x,x,sigma,0,0,0,var_flag1,c) .* A_1st_deriv(x,x,sigma,0,0,0,var_flag2,c))) ./ (A(x,x,sigma,0,0,0, c) .^ 2);        % (x_1,x_2,sigma,u,t,s, var_flag, c)
    F = q_c .* (A_2nd_deriv(x, x, sigma, 0, 0, 0, var_flag1, var_flag2, c) - ...          %(x_1,x_2,sigma,u,t,s, var_flag1, var_flag2, c)
        A_1st_deriv(x,x,sigma,0,0,0,var_flag1,c) .* A_1st_deriv(x,x,sigma,0,0,0,var_flag2,c));        % (x_1,x_2,sigma,u,t,s, var_flag, c)
    % Here we used the fact that A====1 

end

% Now special simple functions for y=z=0
