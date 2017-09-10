% This gives the function which appears in the exponent of the saddle point integral. The represantation here is 
% as a function of s and t rather than of y and z
% we return both the value of the function, and the Jacobian 2nd derivative
% matrix J. J can be 3*3 or 4*4, depending on if we have the numerator or
% denominator.
function [F, J] = FInExponent( x, t, s, sigma, alpha, C_alpha, f_star, JSize, soft_constrain_flag)

TOL = 0.00000000000000001;

TWO = 2; % Should be two !!!

I = complex(0,1); % define i

% Take the principal log value
minus_iy = log(t+1); % -iy
minus_iz = log(s+1)-minus_iy;  % -iz
BIG_J =10;

%%%%F =  -BIG_J * minus_iy^2 + 2*BIG_J*minus_iy * (1-alpha) -BIG_J * minus_iz^2 + 2*BIG_J*minus_iz * alpha * (1-f_star) - ...  % Here its a 'real' 2, not the defined TWO
F =   -TWO .* ZukQuadl(C_alpha, 99, TOL, x, sigma, s, 6) - ...  % g(..)
    TWO .* ZukQuadl(0, C_alpha, TOL, x, sigma, t, 6);       % g(..)

if(soft_constrain_flag)
    F = F  -BIG_J * minus_iy^2 + 2*BIG_J*minus_iy * (1-alpha) -BIG_J * minus_iz^2 + 2*BIG_J*minus_iz * alpha * (1-f_star);  % Here its a 'real' 2, not the defined TWO
else
    F =  F + minus_iy * (1-alpha) + minus_iz * alpha * (1-f_star);
end



% Now calculate the Jacobian.  The parameters order is : x y z f. Note : Here the
% derivatives are with respect to y and z, not t and s. 
J = zeros(4);

if(JSize > 2) % No need to compute J if it is not needed, size <= 2
    J(1,1) = -TWO .* ZukQuadl(C_alpha, 999, TOL, x, sigma, s, 4) - ...
        TWO .* ZukQuadl(0, C_alpha, TOL, x, sigma, t, 4);
    J(1,2) = I*TWO .* ZukQuadl(C_alpha, 999, TOL, x, sigma, s, 2) + ...
        I*TWO .* ZukQuadl(0, C_alpha, TOL, x, sigma, t, 2);
    J(1,3) = I*TWO .* ZukQuadl(C_alpha, 999, TOL, x, sigma, s, 2) ;
    J(1,4) = 0;

%    J(2,1) = I*TWO .* ZukQuadl(C_alpha, 999, TOL, x, sigma, t, 7) + ...
%        I*TWO .* ZukQuadl(0, C_alpha, TOL, x, sigma, s, 7);
    J(2,2) = TWO .* ZukQuadl(C_alpha, 999, TOL, x, sigma, s, 3) + ...
        TWO .* ZukQuadl(0, C_alpha, TOL, x, sigma, t, 3);
    J(2,3) = TWO .* ZukQuadl(C_alpha, 999, TOL, x, sigma, s, 3);
    J(2,4) = 0;

%    J(3,1) = I*TWO .* ZukQuadl(C_alpha, 999, TOL, x, sigma, t, 7);
%    J(3,2) = J(2,3);
    J(3,3) = J(2,3);
    J(3,4) = I*alpha;
    
%    J(4,1) = 0;
%    J(4,2) = 0;
%    J(4,3) = I*alpha;
    J(4,4) = 0;

    % J must be symmetric : 
    J(2,1) = J(1,2);
    J(3,1) = J(1,3);
    J(3,2) = J(2,3);
    J(4,1) = J(1,4);
    J(4,2) = J(2,4);
    J(4,3) = J(3,4);
    
    
    J = J(1:JSize,1:JSize);  % adjust dimension
end
