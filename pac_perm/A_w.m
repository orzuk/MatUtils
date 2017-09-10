% Compute A_w partial derivative of the saddle point integral
function A_w=A_w(x_1,x_2,c,sigma,u,t,s)

A_w=-I*[PInt(x_1,c,sigma).*PInt(x_2,c,sigma).*( u-t-s-1)...
        +PInt(x_1,c,sigma).*(t+1)+PInt(x_2,c,sigma).*(s+1)];
