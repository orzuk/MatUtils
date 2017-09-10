% This gives the function which appears in the exponent of the saddle point integral. The represantation here is 
% as a function of s and t rather than of y and z
% we return both the value of the function, and the Jacobian 2nd derivative
% matrix J. J can be 3*3 or 4*4, depending on if we have the numerator or
% denominator.
function [F, J] = FInExponentTwoSampled( x1, x2, u, t, s, sigma, alpha, C_alpha, f_star, JSize, soft_constrain_flag)

TOL = 0.0000000000000001;

TWO = 2; % Should be two !!!

I = complex(0,1); % define i

% Take the principal log value
minus_iy = log(u+1)-log(s+1); % -iy
minus_iz = log(u+1)-log(t+1);  % -iz
minus_iw = log(t+1)+log(s+1)-log(u+1); % -iw

BIG_J =10;

%%%%F =  -BIG_J * minus_iy^2 + 2*BIG_J*minus_iy * (1-alpha) -BIG_J * minus_iz^2 + 2*BIG_J*minus_iz * alpha * (1-f_star) - ...  % Here its a 'real' 2, not the defined TWO
F =   -TWO .* ZukQuadlManyVars(0, 999, TOL, x1, x2, sigma, u, t, s, 0); % q * log(A)

if(soft_constrain_flag) % Currently not used ...
    F = F  -BIG_J * minus_iy^2 + 2*BIG_J*minus_iy * (1-alpha) -BIG_J * minus_iz^2 + 2*BIG_J*minus_iz * (1-alpha) - ...
        BIG_J * minus_iw^2 + 2*BIG_J*minus_iw * (1-alpha*f_star);  % Here its a 'real' 2, not the defined TWO
else
    F =  F + minus_iy * (1-alpha) + minus_iz * (1-alpha) + minus_iw * (1-alpha*f_star);
end



% Now calculate the Jacobian.  The parameters order is : x y z f. Note : Here the
% derivatives are with respect to y and z, not t and s. 
J = zeros(6);

if(JSize > 2) % No need to compute J if it is not needed, size <= 4
    
%     NOW_PLOT_ALL_JS = 9999
%     x1
%     x2
%     u
%     t
%     s
    FINexponent_sigma = sigma;
    
    % Make all non-zero J
    for i=1:5
        for j=i:5
%            I_is = i
%            J_is = j
            J(i,j) = -2*ZukQuadlManyVars(0, 999, TOL, x1, x2, sigma, u, t, s, 900+10*i+j); % Here take integral from 0 to infinity
%            OOOOF = 12345
        end
    end
    % Make J symmetric
    for i=1:5
        for j=i+1:5
            J(j,i) = J(i,j);
        end
    end
    J(5,6) = alpha*I; J(6,5) = alpha*I; % Here we have F_{fw}

    J = J(1:JSize,1:JSize);  % adjust dimension
end
