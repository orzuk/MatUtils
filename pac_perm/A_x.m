% The partial first derivative of A. One can choose by the flag var_flag
% the variable according to which we differentiate


function    F=A_1st_deriv(x_1,x_2,sigma,u,t,s, var_flag, c)
x_1  = 1; x_2 = 1; y = 3; z = 4;  w = 5; f=6; 

i = sqrt(-1);

switch var_flag 
    case x,
    return P(x_1,c)*P(x_2,c)*(u-t-s) + P_x(x_1,c)*t;
    case y
        return P(x_1,c)*P(x_2,c)*(u-t-s) + P_x(x_1,c)*t;
    case z
        return (-i)*P(x_1,c)*(P(x_2,c)*(u-t)+t+1);
    
end



% Liat's version : 
% function    A_x=A_x(x_1,x_2,c,sigma,u,t,s)
%
% A_x=PInt(x_1,c,sigma) .*PInt(x_2,c,sigma).*( u-t-s)+ ...
%            Px(x_1,c,sigma).*t;
