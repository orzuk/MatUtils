% This shows the weird function that we should integrate in the saddle point equations. 
% The derive flag controls between two kinds of densities, one with P and
% one with P_x
function F = SaddleHelperGaussian_two( c, x,  sigma,u, t,s, deriv_flag)

q_c =  (1./sqrt(2.*pi)) .* exp(-c.^2./2); % Set the density q to be gaussian
A=PInt(x_1,c,sigma).*PInt(x_2,c,sigma)*(u-t-s)...
+PInt(x_1,c,sigma).*t+PInt(x_2,c,sigma).*s+1 ;

% only two-sided flag so far. Our density q is standard gaussian
if(deriv_flag == 1)   % Here take Px
    F =  q_c.*A_x(x_1,x_2,c,sigma,u,t,s)./A;
end

% Here we start with 2nd derivative stuff - The Jacobian : 
if(deriv_flag == 2) % here take Px
    F =  q_c.*A_x(x_2,x_1,c,sigma,u,s,t)./A;
end
if(deriv_flag == 3) % here take P*(1-P)
    F = q_c.*A_y(x_1,x_2,c,sigma,u,t)./A;
end
if(deriv_flag == 4) % here take complicated !!! 
    F = q_c.*A_y(x_2,x_1,c,sigma,u,s)./A;
end
if(deriv_flag == 5) % here take Px without t+1
    F = q_c.*A_w(x_1,x_2,c,sigma,u,t,s)./A;
end

if(deriv_flag == 11) % here take Px without t+1
    A_11=Pxx(x1,c,sigma).*[PInt(x_2,c,sigma).(u-t-s)+t];
    F = (q_c.*A_11.*A-A_x(x_1,x_2,c,sigma, u,t,s).^2)/A;
end


% This is a special flag for computing g. Log and no denominator ! 

%     A_11=Pxx(x1,c,sigma).*[PInt(x_2,c,sigma).(u-t-s)+t];
