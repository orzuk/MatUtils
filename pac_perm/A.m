% The function A, which is used to compute the std. of f
function    F=A(x_1,x_2,sigma,u,t,s, c)

% x_1
% c
% sigma
% F1 = PInt(x_1,c,sigma)
% F2 = PInt(x_2,c,sigma)
% F3 = (u-t-s)
% F4 = PInt(x_1,c,sigma)
% F5 = t
% F6 = PInt(x_2,c,sigma)
% F7 = 1
% 
% input1 = (x_1-c) / sigma
% input2 = (x_2-c) / sigma
% 
F = PInt(x_1,c,sigma).*PInt(x_2,c,sigma) .* (u-t-s) + PInt(x_1,c,sigma).*t + PInt(x_2,c,sigma).*s + 1;

