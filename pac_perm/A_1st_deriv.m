% The partial first derivative of A. One can choose by the flag var_flag
% the variable according to which we differentiate


function    F=A_1st_deriv(x_1,x_2,sigma,u,t,s, var_flag, c)
%x_1 = 1; x_2 = 2; y = 3; z = 4;  w = 5; f = 6; 


I = sqrt(-1);

switch var_flag 
    case 1 % x_1
        F =  Px(x_1,c,sigma).*( PInt(x_2,c,sigma).*(u-t-s) + t );   % New - Corrected
    case 2 % x_2
        F =  Px(x_2,c,sigma).*( PInt(x_1,c,sigma).*(u-t-s) + s );  % New - Corrected
    case 3 % y
        F =  (-I).*PInt(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t)+t+1);
    case 4 % z
        F =  (-I).*PInt(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-s)+s+1);
    case 5 % w
        F =  (-I).*[PInt(x_1,c,sigma).*PInt(x_2,c,sigma).*(u-t-s-1) + PInt(x_1,c,sigma).*(t+1) + PInt(x_2,c,sigma).*(s+1)];
end

% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % Try to already put zeros
% % % % % switch var_flag 
% % % % %     case x_1,
% % % % %         F =  0;
% % % % %     case x_2
% % % % %         F =  0; 
% % % % %     case y
% % % % %         F =  (-I).*PInt(x_1,c,sigma);
% % % % %     case z
% % % % %         F =  (-I).*PInt(x_2,c,sigma);
% % % % %     case w
% % % % %         F =  (-I).*[-PInt(x_1,c,sigma).*PInt(x_2,c,sigma) + PInt(x_1,c,sigma) + PInt(x_2,c,sigma)];
% % % % % end
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Liat's version : 
% function    A_x=A_x(x_1,x_2,c,sigma,u,t,s)
%
% A_x=PInt(x_1,c,sigma) .*PInt(x_2,c,sigma).*( u-t-s)+ ...
%            Px(x_1,c,sigma).*t;
