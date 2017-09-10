% A lot of calculations for computing the std, in the gaussian case for now ..
function inf_limit_std = compute_inf_limit_std(rand_flag, one_side_flag, true_corr_flag, sigma, alpha, C_alpha, x_alpha, f_star)
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
TOL = 0.00000000001;
I = sqrt(-1);

TWO =2;

if(true_corr_flag == TRUE_AND_SAMPLED)
    % New ! (28.6.05) - Compute also the (approx.) standard
    % deviation analytically
    a = zeros(4);
    a(2) = 2 * I * SaddleHelperGaussianQuadl(C_alpha, 999, TOL, x_alpha, sigma, 0, 2);
    a(1) = 2 * I * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_alpha, sigma, 0, 2) + a(2);

    a(4) = 2 * SaddleHelperGaussianQuadl(C_alpha, 999, TOL, x_alpha, sigma, 0, 3);
    a(3) = 2 * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_alpha, sigma, 0, 3) + a(4);

    inf_limit_std = sqrt( a(4) *(a(1)-a(2))^2 + (a(3) - a(4)) * a(2)^2 ) / (alpha * a(1));



% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%%% Computation of J is not really needed - it is just for debugging
% % % % %     % Compute J at the saddle-point
% % % % %     J = zeros(4);
% % % % % 
% % % % %     J(1,1) = -TWO .* SaddleHelperGaussianQuadl(C_alpha, 999, TOL, x_alpha, sigma, 0, 4) - ...
% % % % %         TWO .* SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_alpha, sigma, 0, 4);
% % % % %     J(1,2) = I*TWO .* SaddleHelperGaussianQuadl(C_alpha, 999, TOL, x_alpha, sigma, 0, 2) + ...
% % % % %         I*TWO .* SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_alpha, sigma, 0, 2);
% % % % %     J(1,3) = I*TWO .* SaddleHelperGaussianQuadl(C_alpha, 999, TOL, x_alpha, sigma, 0, 2) ;
% % % % %     J(1,4) = 0;
% % % % % 
% % % % %     J(2,2) = TWO .* SaddleHelperGaussianQuadl(C_alpha, 999, TOL, x_alpha, sigma, 0, 3) + ...
% % % % %         TWO .* SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_alpha, sigma, 0, 3);
% % % % %     J(2,3) = TWO .* SaddleHelperGaussianQuadl(C_alpha, 999, TOL, x_alpha, sigma, 0, 3);
% % % % %     J(2,4) = 0;
% % % % % 
% % % % %     J(3,3) = J(2,3);
% % % % %     J(3,4) = I*alpha;
% % % % %     J(4,4) = 0;
% % % % % 
% % % % %     % J must be symmetric :
% % % % %     J(2,1) = J(1,2);
% % % % %     J(3,1) = J(1,3);
% % % % %     J(3,2) = J(2,3);
% % % % %     J(4,1) = J(1,4);
% % % % %     J(4,2) = J(2,4);
% % % % %     J(4,3) = J(3,4);
% % % % % 
% % % % %     J;
% % % % %     TRUE_AND_SAMP_INV_J_IS = inv(J);
% % % % %     TRUE_AND_SAMP_INV_J3_IS = inv(J(1:3,1:3));
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%     inf_limit_std3 = sqrt(det(J) / (det(J(1:3,1:3))));
% 
%     inf_limit_std4 = 1/inf_limit_std3;
% 
%     inf_limit_std5 = sqrt(1) / (alpha*sqrt(TRUE_AND_SAMP_INV_J3_IS(3,3)));
% 
%     %%% TRUE_AND_SAMP_second_diff_should_be_zero = inf_limit_std-inf_limit_std3 %%%    This one is bad
%     TRUE_AND_SAMP_third_diff_should_be_zero = inf_limit_std-inf_limit_std4
%     TRUE_AND_SAMP_FOURTH_diff_should_be_zero = inf_limit_std-inf_limit_std5


else % here the major calculations, between two sampled distributions

    J = zeros(6);

    % Make all non-zero J
    for i=1:5
        for j=i:5
            %%%            J(i,j) = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*i+j); % Here take integral from 0 to infinity
            J(i,j) = -2*SaddleHelperGaussianQuadlManyVars(0, 999, TOL, x_alpha, x_alpha, sigma, 0, 0, 0, 900+10*i+j);
        end
    end
    % Make J symmetric
    for i=1:5
        for j=i+1:5
            J(j,i) = J(i,j);
        end
    end
    J(5,6) = alpha*I; J(6,5) = alpha*I; % Here we have F_{fw}





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  Numerical Computations for Debugging :
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % %     % Now calculate J numerically!
    % % % % %     NUM_J = zeros(6);
    % % % % %
    % % % % %
    % % % % %     % Calculate the 2-nd derivatives numerically:
    % % % % %     % Start by the diagonal - F_xx
    % % % % %     delta_small = 0.0001; % The delta we use for calculating the 2nd derivatives numerically
    % % % % %     saddle_v_vec = [x_alpha,x_alpha,0,0,0,f_star];
    % % % % %     for i=1:6
    % % % % %
    % % % % %         delta_small_vec = zeros(6,1); delta_small_vec(i) = delta_small;
    % % % % %         % We need to extract t,s,u from y,z,w :
    % % % % %         % This is the opposite:
    % % % % %         %         minus_iy = log(u+1)-log(s+1); % -iy
    % % % % %         %         minus_iz = log(u+1)-log(t+1);  % -iz
    % % % % %         %         minus_iw = log(t+1)+log(s+1)-log(u+1); % -iw
    % % % % %         % This is the correct way:
    % % % % %         t_plus = exp(-I*(delta_small_vec(3)+delta_small_vec(5)))-1;
    % % % % %         s_plus = exp(-I*(delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %         u_plus = exp(-I*(delta_small_vec(3)+delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %
    % % % % %         delta_small_vec(i) = -delta_small;
    % % % % %         t_minus = exp(-I*(delta_small_vec(3)+delta_small_vec(5)))-1;
    % % % % %         s_minus = exp(-I*(delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %         u_minus = exp(-I*(delta_small_vec(3)+delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %
    % % % % %         % Now change the vec
    % % % % %         %delta_small_vec = zeros(6,1)+delta_small;
    % % % % %         %delta_small_vec(3) = u_is; delta_small_vec(4) = t_is; delta_small_vec(5) = s_is;
    % % % % %
    % % % % %         % FInExponentTwoSampled( x1, x2, u, t, s, sigma, alpha, C_alpha, f_star, JSize, soft_constrain_flag)
    % % % % %         delta_plus_vec = saddle_v_vec; delta_plus_vec(i) = delta_plus_vec(i) + delta_small;
    % % % % %         delta_minus_vec = saddle_v_vec; delta_minus_vec(i) = delta_plus_vec(i) - delta_small;
    % % % % %
    % % % % %         [A_F DUM_J] = FInExponentTwoSampled(delta_plus_vec(1), delta_plus_vec(2), u_plus, t_plus, ...
    % % % % %             s_plus, sigma, alpha, C_alpha, delta_plus_vec(6), 1, 0);
    % % % % %         [B_F DUM_J] = FInExponentTwoSampled( delta_minus_vec(1), delta_minus_vec(2), u_minus, t_minus, ...
    % % % % %             s_minus, sigma, alpha, C_alpha, delta_minus_vec(6), 1, 0);
    % % % % %         [C_F DUM_J] = FInExponentTwoSampled( x_alpha, x_alpha, 0, 0, 0, sigma, alpha, C_alpha, f_star, 1, 0);
    % % % % %
    % % % % %         NUM_J(i,i) = ( A_F + B_F - 2*C_F ) ./ (delta_small .^2)
    % % % % %
    % % % % %     end
    % % % % %
    % % % % %     % Now work on the non-diagonal elements :
    % % % % %     for i=1:6
    % % % % %         for j=i+1:6
    % % % % %             delta_small_vec = zeros(6,1); delta_small_vec(i) = delta_small; delta_small_vec(j) = delta_small;
    % % % % %             % We need to extract t,s,u from y,z,w :
    % % % % %             % This is the opposite:
    % % % % %             %         minus_iy = log(u+1)-log(s+1); % -iy
    % % % % %             %         minus_iz = log(u+1)-log(t+1);  % -iz
    % % % % %             %         minus_iw = log(t+1)+log(s+1)-log(u+1); % -iw
    % % % % %             % This is the correct way:
    % % % % %
    % % % % %
    % % % % %             t_plus = exp(-I*(delta_small_vec(3)+delta_small_vec(5)))-1;
    % % % % %             s_plus = exp(-I*(delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %             u_plus = exp(-I*(delta_small_vec(3)+delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %
    % % % % %             delta_small_vec(i) = -delta_small; delta_small_vec(j) = -delta_small;
    % % % % %             t_minus = exp(-I*(delta_small_vec(3)+delta_small_vec(5)))-1;
    % % % % %             s_minus = exp(-I*(delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %             u_minus = exp(-I*(delta_small_vec(3)+delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %
    % % % % %
    % % % % %             delta_plus_vec = saddle_v_vec; delta_plus_vec(i) = delta_plus_vec(i) + delta_small; delta_plus_vec(j) = delta_plus_vec(j) + delta_small;
    % % % % %             delta_minus_vec = saddle_v_vec; delta_minus_vec(i) = delta_minus_vec(i) - delta_small; delta_minus_vec(j) = delta_minus_vec(j) - delta_small;
    % % % % %
    % % % % %
    % % % % %             [A_F DUM_J] = FInExponentTwoSampled(delta_plus_vec(1), delta_plus_vec(2), u_plus, t_plus, ...
    % % % % %                 s_plus, sigma, alpha, C_alpha, delta_plus_vec(6), 1, 0);
    % % % % %             [B_F DUM_J] = FInExponentTwoSampled( delta_minus_vec(1), delta_minus_vec(2), u_minus, t_minus, ...
    % % % % %                 s_minus, sigma, alpha, C_alpha, delta_minus_vec(6), 1, 0);
    % % % % %             [C_F DUM_J] = FInExponentTwoSampled( x_alpha, x_alpha, 0, 0, 0, sigma, alpha, C_alpha, f_star, 1, 0);
    % % % % %
    % % % % %
    % % % % %             NUM_J(i,j) = ( A_F + B_F - 2*C_F ) ./ (2 .* delta_small.^2) - NUM_J(i,i) ./ 2 - NUM_J(j,j) ./ 2;
    % % % % %
    % % % % %         end
    % % % % %     end
    % % % % %
    % % % % %
    % % % % %      for i=1:6
    % % % % %         for j=i+1:6
    % % % % %             NUM_J(j,i) = NUM_J(i,j);
    % % % % %         end
    % % % % %      end

    % Now compare :


    %%%J = NUM_J; % Wild try
    % After calculating J, we need to invert it :
    i_alpha_vec = zeros(5, 1); i_alpha_vec(end) = -I*alpha;  % We assume delta F is one here !!
    delta_v = inv(J(1:5,1:5)) * i_alpha_vec;
    delta_v(end + 1) = 1; % Put delta F as one !

    INV_J_IS = inv(J);
    INV_J5_IS = inv(J(1:5,1:5));

    inf_limit_std2 = 1.0/sqrt(transpose(delta_v) * J * delta_v); % See if this is complex or not !


    % Here we already directly inverted J (symbolically), and only need to
    % calculate the relevant F_xy derivatives
    % x1 = 1 x2 = 2 y = 3 z = 4 w = 5 f = 6
    F_x1y = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*1+3);
    F_x1w = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*1+5);
    F_x2z = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*2+4);
    F_x2w = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*2+5);
    F_yy = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*3+3);
    F_yw = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*3+5);
    F_zz = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*4+4);
    F_zw = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*4+5);
    F_ww = -2*SaddleHelperGaussianQuadl(0, 999, TOL, x_alpha, sigma, 0, 900+10*5+5);


    inf_limit_variance = 1.0*(F_ww*F_x1y^2*F_x2z^2-2*F_x1w*F_yw*F_x1y*F_x2z^2-2*F_x2w*F_zw*F_x1y^2*F_x2z+F_x1w^2*F_x2z^2*F_yy+F_x2w^2*F_x1y^2*F_zz)...
        /(alpha^2*F_x1y^2*F_x2z^2);   % Should be real > 0 !!!!!
    inf_limit_std = sqrt(1.0*(F_ww*F_x1y^2*F_x2z^2-2*F_x1w*F_yw*F_x1y*F_x2z^2-2*F_x2w*F_zw*F_x1y^2*F_x2z+F_x1w^2*F_x2z^2*F_yy+F_x2w^2*F_x1y^2*F_zz))/...
        abs(alpha*F_x1y*F_x2z); % Try only denominator - This must be the square-root of the previous !!!!


    % Now try a third approach. Simply divide determinantes :
%     inf_limit_std3 = sqrt(det(J) / (det(J(1:5,1:5))));
% 
%     inf_limit_std4 = 1/inf_limit_std3;
% 
%     inf_limit_std5 = sqrt(1.0) / (alpha*sqrt(INV_J5_IS(5,5)));
% 
%     % Check if we inverted all matrices correctly ...
% %     TWO_SAMPS_diff_should_be_zero = inf_limit_std-inf_limit_std2
% %     %%    TWO_SAMPS_second_diff_should_be_zero = inf_limit_std-inf_limit_std3 %%%    This one is bad
% %     TWO_SAMPS_third_diff_should_be_zero = inf_limit_std-inf_limit_std4
% %     TWO_SAMPS_FOURTH_diff_should_be_zero = inf_limit_std-inf_limit_std5

    %%%inf_limit_std = inf_limit_std4;

    % Maple Computed  the following :
    %sigma :=
    %sqrt(F_ww*F_x1y^2*F_x2z^2-2*F_x1w*F_yw*F_x1y*F_x2z^2-2*F_x2w*F_zw*F_x1y^2*F_x2z+F_x1w^2*F_x2z^2*F_yy+F_x2w^2*F_x1y^2*F_zz)* .
    % sqrt(2)*sqrt(1/(alpha^2*F_x1y^2*F_x2z^2*N))




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Gradient Numerical Computations for debugging:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % %     % Now calculate the gradient. Maybe this is good ? The grad is supposed
    % % % % %     % to be zero in both calculations !!
    % % % % %     G = zeros(6,1); NUM_G = zeros(6,1);
    % % % % %
    % % % % %     for i=1:5
    % % % % %
    % % % % %         G(i) = -2*SaddleHelperGaussianQuadlManyVars(0, 999, TOL, x_alpha, x_alpha, sigma, 0, 0, 0, 80+i);
    % % % % %
    % % % % %         GGG_i = G(i)
    % % % % %         if(i==3)
    % % % % %             G(i) = G(i) - I*(1-alpha);
    % % % % %         end
    % % % % %         if(i==4)
    % % % % %             G(i) = G(i) - I*(1-alpha);
    % % % % %         end
    % % % % %         if(i==5)
    % % % % %             G(i) = G(i) - I*(1-alpha*f_star);
    % % % % %         end
    % % % % %     end
    % % % % %     G(6) = alpha*I*0;
    % % % % %
    % % % % %     for i=1:6
    % % % % %           delta_small_vec = zeros(6,1); delta_small_vec(i) = delta_small;
    % % % % %         % We need to extract t,s,u from y,z,w :
    % % % % %         % This is the opposite:
    % % % % %         %         minus_iy = log(u+1)-log(s+1); % -iy
    % % % % %         %         minus_iz = log(u+1)-log(t+1);  % -iz
    % % % % %         %         minus_iw = log(t+1)+log(s+1)-log(u+1); % -iw
    % % % % %         % This is the correct way:
    % % % % %         t_plus = exp(-I*(delta_small_vec(3)+delta_small_vec(5)))-1;
    % % % % %         s_plus = exp(-I*(delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %         u_plus = exp(-I*(delta_small_vec(3)+delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %
    % % % % %         delta_small_vec(i) = -delta_small;
    % % % % %         t_minus = exp(-I*(delta_small_vec(3)+delta_small_vec(5)))-1;
    % % % % %         s_minus = exp(-I*(delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %         u_minus = exp(-I*(delta_small_vec(3)+delta_small_vec(4)+delta_small_vec(5)))-1;
    % % % % %
    % % % % %         % Now change the vec
    % % % % %         %delta_small_vec = zeros(6,1)+delta_small;
    % % % % %         %delta_small_vec(3) = u_is; delta_small_vec(4) = t_is; delta_small_vec(5) = s_is;
    % % % % %
    % % % % %         % FInExponentTwoSampled( x1, x2, u, t, s, sigma, alpha, C_alpha, f_star, JSize, soft_constrain_flag)
    % % % % %         delta_plus_vec = saddle_v_vec; delta_plus_vec(i) = delta_plus_vec(i) + delta_small;
    % % % % %         delta_minus_vec = saddle_v_vec; delta_minus_vec(i) = delta_plus_vec(i) - delta_small;
    % % % % %
    % % % % %         [A_F DUM_J] = FInExponentTwoSampled(delta_plus_vec(1), delta_plus_vec(2), u_plus, t_plus, ...
    % % % % %             s_plus, sigma, alpha, C_alpha, delta_plus_vec(6), 1, 0);
    % % % % %         [B_F DUM_J] = FInExponentTwoSampled( delta_minus_vec(1), delta_minus_vec(2), u_minus, t_minus, ...
    % % % % %             s_minus, sigma, alpha, C_alpha, delta_minus_vec(6), 1, 0);
    % % % % %
    % % % % %         NUM_G(i) = ( B_F - A_F ) ./ (2.*delta_small);
    % % % % %
    % % % % %     end
    % % % % %
    % % % % %
    % % % % %     G
    % % % % %     NUM_G
    % % % % %     GRAD_DIFF = G-NUM_G
    % % % % %
    % % % % % end


end



% Now perform a correction to avoid too small and imaginary st.d.'s
if(abs(inf_limit_std) < TOL)
    inf_limit_std = TOL;
end



