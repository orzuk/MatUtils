% Compute A_y partial derivative of the saddle point integral
function A_y=A_y(x_1,x_2,c,sigma,u,t)

A_y=-I*PInt(x_1,c,sigma).*[PInt(x_2,c,sigma).*( u-t)+t+1];
