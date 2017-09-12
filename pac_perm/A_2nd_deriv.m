% The partial second derivative of A. One can choose by the two flags
% var_flag1, var_flag2 the variables according to which we differentiate
function    F=A_2nd_deriv(x_1,x_2,sigma,u,t,s, var_flag1, var_flag2, c)
y = 3; z = 4;  w = 5; f = 6; 


I = sqrt(-1);


switch var_flag1 
    case 1, % x_1
        switch var_flag2
            case 1, % x_1
                F =  Pxx(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t-s) + t);
            case 2, % x_2
                F =  Px(x_1,c,sigma).*Px(x_2,c,sigma).*(u-t-s);
            case y,
                F =  (-I).*Px(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t)+t+1);
            case z,
                F =  (-I).*Px(x_1,c,sigma).*PInt(x_2,c,sigma).*(u-s);
            case w,
                F =  (-I).*Px(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t-s-1)+t+1); 
        end
    case 2 % x_2
        switch var_flag2
            case 2, % x_2
                F =  Pxx(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-t-s) + s);
            case y,
                F =  (-I).*Px(x_2,c,sigma).*PInt(x_1,c,sigma).*(u-t);
            case z,
                F =  (-I).*Px(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-s)+s+1);
            case w,
                F =  (-I).*Px(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-t-s-1)+s+1);
        end
    case y
        switch var_flag2
            case y,
                F =  -PInt(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t)+t+1);
            case z,
                F =  -PInt(x_1,c,sigma).*PInt(x_2,c,sigma).*(u+1);  % Here need to check if it is u+1 or u-1 !!! Should be '+' 
            case w,
                F =  -PInt(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t)+t+1);
        end
    case z   % here A_zz = A_zw - Still ??? 
        F =  -PInt(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-s)+s+1);
    case w % here we are left with only A_ww
        F =  -PInt(x_1,c,sigma).*PInt(x_2,c,sigma).*(u-t-s-1) - PInt(x_1,c,sigma).*(t+1) - PInt(x_2,c,sigma).*(s+1);
end


% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % Try to already put zeros
% % % % % switch var_flag1 
% % % % %     case x_1,
% % % % %         switch var_flag2
% % % % %             case x_1,
% % % % %                 F =  0;
% % % % %             case x_2,
% % % % %                 F =  0;
% % % % %             case y,
% % % % %                 F =  (-i).*Px(x_1,c,sigma);
% % % % %             case z,
% % % % %                 F =  0;
% % % % %             case w,
% % % % %                 F =  (-i).*Px(x_1,c,sigma).*(1-PInt(x_2,c,sigma)); 
% % % % %         end
% % % % %     case x_2
% % % % %         switch var_flag2
% % % % %             case x_2,
% % % % %                 F = 0;
% % % % %             case y,
% % % % %                 F = 0;
% % % % %             case z,
% % % % %                 F =  (-i).*Px(x_2,c,sigma);
% % % % %             case w,
% % % % %                 F =  (-i).*Px(x_2,c,sigma).*(1-PInt(x_1,c,sigma));
% % % % %         end
% % % % %     case y
% % % % %         switch var_flag2
% % % % %             case y,
% % % % %                 F =  -PInt(x_1,c,sigma);
% % % % %             case z,
% % % % %                 F =  -PInt(x_1,c,sigma).*PInt(x_2,c,sigma);  % Here need to check if it is u+1 or u-1 !!! Should be '+' 
% % % % %             case w,
% % % % %                 F =  -PInt(x_1,c,sigma);
% % % % %         end
% % % % %     case z   % here A_zz = A_zw - Still ??? 
% % % % %         F =  -PInt(x_2,c,sigma);
% % % % %     case w % here we are left with only A_ww
% % % % %         F =  PInt(x_1,c,sigma).*PInt(x_2,c,sigma) - PInt(x_1,c,sigma) - PInt(x_2,c,sigma);
% % % % % end
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



