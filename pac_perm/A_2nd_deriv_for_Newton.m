% The partial second derivative of A. One can choose by the two flags
% var_flag1, var_flag2 the variables according to which we differentiate
function    F=A_2nd_deriv_for_Newton(x_1,x_2,sigma,u,t,s, var_flag1, var_flag2, c)
y = 3; z = 4;  w = 5; f = 6; % y is like u, z is t, w is like s ...


I = sqrt(-1);


switch var_flag1
    case x_1,
        switch var_flag2
            case x_1,
                F =  Pxx(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t-s) + t);
            case x_2,
                F =  Px(x_1,c,sigma).*Px(x_2,c,sigma).*(u-t-s);
            case y, % u
                F =  Px(x_1,c,sigma).* PInt(x_2,c,sigma);
            case z, % t
                F =  Px(x_1,c,sigma).*(1-PInt(x_2,c,sigma));
            case w, % s
                F =  -Px(x_1,c,sigma).* PInt(x_2,c,sigma);
        end
    case x_2
        switch var_flag2
            case x_1,
                F =  Px(x_1,c,sigma).*Px(x_2,c,sigma).*(u-t-s); 
            case x_2,
                F =  Pxx(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-t-s) + s);
            case y, % u
                F =  Px(x_2,c,sigma).*PInt(x_1,c,sigma);
            case z, % t
                F =  -Px(x_2,c,sigma).*PInt(x_1,c,sigma);
            case w, % s
                F =  Px(x_2,c,sigma).*(1-PInt(x_1,c,sigma));
        end
    case y
        switch var_flag2
            case x_1, 
                F =  (-I).*Px(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t)+t+1);
            case x_2, 
                F =  (-I).*Px(x_2,c,sigma).*PInt(x_1,c,sigma).*(u-t);
            case y,
                F =  (-I).*PInt(x_1,c,sigma).*PInt(x_2,c,sigma);
            case z,
                F =  I.*PInt(x_1,c,sigma).*(PInt(x_2,c,sigma)-1);  % Here need to check if it is u+1 or u-1 !!! Should be '+'
            case w,
                F = 0;
        end
    case z   % here A_zz = A_zw - Still ???
        switch var_flag2
            case x_1, 
                F =  (-I).*Px(x_1,c,sigma).*PInt(x_2,c,sigma).*(u-s);
            case x_2
                F =  (-I).*Px(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-s)+s+1);
            case y,
                F = -I.*PInt(x_1,c,sigma).*PInt(x_2,c,sigma);
            case z,
                F = 0;
            case w,
                F = I.*PInt(x_2,c,sigma).*(PInt(x_1,c,sigma)-1);
        end
    case w % here we are left with only A_ww
        switch var_flag2
            case x_1,
                F =  (-I).*Px(x_1,c,sigma).*(PInt(x_2,c,sigma).*(u-t-s-1)+t+1);
            case x_2
                F =  (-I).*Px(x_2,c,sigma).*(PInt(x_1,c,sigma).*(u-t-s-1)+s+1);
            case y,
                F = -I.*PInt(x_1,c,sigma).*PInt(x_2,c,sigma);
            case z,
                F = I.*PInt(x_1,c,sigma).* (PInt(x_2,c,sigma)-1);
            case w,
                F = I.*PInt(x_2,c,sigma).* (PInt(x_1,c,sigma)-1);
        end
end





