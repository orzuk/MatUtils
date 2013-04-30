% The partial first derivative of A. One can choose by the flag var_flag
% the variable according to which we differentiate
% Here we differentiate with respect to x1,x2,u,t,s

function    F=A_1st_deriv_for_Newton(x_1,x_2,sigma,u,t,s, var_flag, c)
y = 3; z = 4;  w = 5; f = 6; 


I = sqrt(-1);

switch var_flag 
    case 1, % x_1
        F =  Px(x_1,c,sigma).*( PInt(x_2,c,sigma).*(u-t-s) + t );   % New - Corrected
    case 2, % x_2
        F =  Px(x_2,c,sigma).*( PInt(x_1,c,sigma).*(u-t-s) + s );  % New - Corrected
    case y % u
        F =  PInt(x_1,c,sigma).* PInt(x_2,c,sigma);
    case z % t
        F =  PInt(x_1,c,sigma).* (1-PInt(x_2,c,sigma));
    case w % s
        F =  PInt(x_2,c,sigma).* (1-PInt(x_1,c,sigma));
end




% Liat's version : 
% function    A_x=A_x(x_1,x_2,c,sigma,u,t,s)
%
% A_x=PInt(x_1,c,sigma) .*PInt(x_2,c,sigma).*( u-t-s)+ ...
%            Px(x_1,c,sigma).*t;
