% Compute first order saddle point approximation of the delta distribution for the overlap f
% This should give a good approximation for for large but finite Ngenes
function [saddle_kept_frac_dist sigma_expansion] = ...
    compute_saddle_point_first_order_estimation(rand_flag, one_side_flag, true_corr_flag, dist_std, ...
    nsamples, num_vars, alpha, f_star, x_alpha , res, soft_constrain_flag)


saddle_kept_frac_dist = zeros(1,1/res-1);
f_vec = [res:res:1-res]; % avoid zero and one to avoid troubles ...

UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
TOL = 0.0000001; % This is just TOL for probabilities ..

I = sqrt(-1);


%Note :  Here we must calculate for every f_star, not just the 'correct'
%one !!!


% Get the C_alpha 
if(one_side_flag)
    if(rand_flag == GAUSSIAN)   
        C_alpha = norminv(1-alpha);  
    end
    if(rand_flag == UNIFORM)   
        C_alpha = sqrt(3) * (1-2*alpha);  
    end
    if(rand_flag == LINEAR)   
        C_alpha = sqrt(6) - 2*sqrt(3*alpha);  
    end
else % Take two sides
    if(rand_flag == GAUSSIAN)   
        C_alpha = norminv(1-0.5*alpha);  
    end
    if(rand_flag == UNIFORM)   
        C_alpha = sqrt(3) * (1-alpha);  
    end
    if(rand_flag == LINEAR)   
        C_alpha = sqrt(6) - 2*sqrt(1.5*alpha);  
    end
end



sigma = 1 / (dist_std*sqrt(nsamples)); % st.d. is proportional to 1/sqrt(N)
SIGMA_IS = sigma



if(true_corr_flag==TRUE_AND_SAMPLED)   % here we compare the true and measured correlations
    if(one_side_flag)
        if(rand_flag == UNIFORM)

        end
        if(rand_flag == LINEAR)

        end
        if(rand_flag == GAUSSIAN)

        end
    else  % here take two sides
        if(rand_flag == UNIFORM)

        end
        if(rand_flag == LINEAR)

        end
        if(rand_flag == GAUSSIAN)
            % Note : For now we do only two-sides with gaussian dists.
            % Initilize the X solution
            x_vec = [x_alpha, 0, 0]';
            %prev_x_vec; % A linear extrapulation trick to speed calculations
            [THIS_F_SHOULD_BE_ZERO J4] = FInExponent( x_vec(1), x_vec(2), x_vec(3), sigma, alpha, C_alpha, f_star, 4,  soft_constrain_flag); 
            
            det_2nd_deriv_should_be_positive = det(J4);
            
            % find the index of the correct f_star
            f_star_start = sum(f_vec < f_star);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          Plot six planes 
           npoints = 1;
           y_loc_vec = [-npoints*res:res:npoints*res]; z_loc_vec = y_loc_vec;
           x_loc_vec = x_alpha + y_loc_vec; f_loc_vec = f_star + y_loc_vec;
           
           
%            figure; hold on; subplot(2,3,1);
%            
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            % Do XY plane
%            loc_F_mat = zeros(length(y_loc_vec));
%            for j=1:length(x_loc_vec)
%                now_j = j
%                for k=1:length(y_loc_vec)   
%                   [temp_F J_temp] = FInExponent( x_loc_vec(j), y_loc_vec(k), 0, sigma, alpha, C_alpha, f_star, 2);
%                   loc_F_mat(j, k) = temp_F;
%                end
%            end           
%            subplot(2,3,1); hold on; imagesc(loc_F_mat); colorbar; title('F near saddle in XY plane'); 
%             
%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            % Do XZ plane
%            loc_F_mat = zeros(length(y_loc_vec));
%            for j=1:length(x_loc_vec)
%                now_j = j
%                for k=1:length(z_loc_vec)   
%                   [temp_F J_temp] = FInExponent( x_loc_vec(j), 0, z_loc_vec(k), sigma, alpha, C_alpha, f_star, 2);
%                   loc_F_mat(j, k) = temp_F;
%                end
%            end           
%            subplot(2,3,2); hold on; imagesc(loc_F_mat); colorbar; title('F near saddle in XZ plane'); 
%            
%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            % Do Xf plane
%            loc_F_mat = zeros(length(y_loc_vec));
%            for j=1:length(x_loc_vec)
%                now_j = j
%                for k=1:length(f_loc_vec)   
%                   [temp_F J_temp] = FInExponent( x_loc_vec(j), 0, 0, sigma, alpha, C_alpha, f_loc_vec(k), 2);
%                   loc_F_mat(j, k) = temp_F;
%                end
%            end           
%            subplot(2,3,3); hold on; imagesc(loc_F_mat); colorbar; title('F near saddle in Xf plane'); 
%            
%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            % Do YZ plane
%            loc_F_mat = zeros(length(y_loc_vec));
%            for j=1:length(y_loc_vec)
%                now_j = j
%                for k=1:length(z_loc_vec)   
%                   [temp_F J_temp] = FInExponent(x_alpha, y_loc_vec(j),  z_loc_vec(k), sigma, alpha, C_alpha, f_star, 2);
%                   loc_F_mat(j, k) = temp_F;
%                end
%            end           
%            subplot(2,3,4); hold on; imagesc(loc_F_mat); colorbar; title('F near saddle in YZ plane'); 
%            
%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            % Do Yf plane
%            loc_F_mat = zeros(length(y_loc_vec));
%            for j=1:length(y_loc_vec)
%                now_j = j
%                for k=1:length(y_loc_vec)   
%                   [temp_F J_temp] = FInExponent( x_alpha, y_loc_vec(j), 0, sigma, alpha, C_alpha, f_loc_vec(k), 2);
%                   loc_F_mat(j, k) = temp_F;
%                end
%            end           
%            subplot(2,3,5); hold on; imagesc(loc_F_mat); colorbar; title('F near saddle in Yf plane'); 
%            
%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            % Do Zf plane
%            loc_F_mat = zeros(length(y_loc_vec));
%            for j=1:length(z_loc_vec)
%                now_j = j
%                for k=1:length(f_loc_vec)   
%                   [temp_F J_temp] = FInExponent( x_alpha, 0, z_loc_vec(j), sigma, alpha, C_alpha, f_loc_vec(k), 2);
%                   loc_F_mat(j, k) = temp_F;
%                end
%            end           
%            subplot(2,3,6); hold on; imagesc(loc_F_mat); colorbar; title('F near saddle in Zf plane'); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           % All this was nice but now we want to find if we have here a
           % maximum. For this, we move f_star a little bit, and draw the
           % YZ plane for different values of X : 
           
           figure; hold on; 

           % Do YZ plane
           loc_F_mat = zeros(length(y_loc_vec));
           for j=1:length(y_loc_vec)
               now_j = j;
               for k=1:length(z_loc_vec)   
                  [temp_F J_temp] = FInExponent(x_alpha+res*0.24, y_loc_vec(j),  z_loc_vec(k), sigma, alpha, C_alpha, f_star-0.2*res, 2, soft_constrain_flag);
                  loc_F_mat(j, k) = temp_F;
               end
           end


           % Do YZ plane
           loc_F_mat2 = zeros(length(y_loc_vec));
           for j=1:length(y_loc_vec)
               now_j = j;
               for k=1:length(z_loc_vec)
                   [temp_F J_temp] = FInExponent(x_alpha+res*0.08, y_loc_vec(j),  z_loc_vec(k), sigma, alpha, C_alpha, f_star-0.2*res, 2, soft_constrain_flag);
                   loc_F_mat2(j, k) = temp_F;
               end
           end


           % Do YZ plane
           loc_F_mat3 = zeros(length(y_loc_vec));
           for j=1:length(y_loc_vec)
               now_j = j;
               for k=1:length(z_loc_vec)
                   [temp_F J_temp] = FInExponent(x_alpha-res*0.08, y_loc_vec(j),  z_loc_vec(k), sigma, alpha, C_alpha, f_star-0.2*res, 2, soft_constrain_flag);
                   loc_F_mat3(j, k) = temp_F;
               end
           end


           % Do YZ plane
           loc_F_mat4 = zeros(length(y_loc_vec));
           for j=1:length(y_loc_vec)
               now_j = j;
               for k=1:length(z_loc_vec)
                   [temp_F J_temp] = FInExponent(x_alpha-res*0.24, y_loc_vec(j),  z_loc_vec(k), sigma, alpha, C_alpha, f_star-0.2*res, 2, soft_constrain_flag);
                   loc_F_mat4(j, k) = temp_F;
               end
           end

           low = min(min(loc_F_mat)); low = min(low, min(min(loc_F_mat2))); low = min(low, min(min(loc_F_mat3))); low = min(low, min(min(loc_F_mat4)));
           high = max(max(loc_F_mat)); high = max(low, max(max(loc_F_mat2))); high = max(low, max(max(loc_F_mat3))); high = max(low, max(max(loc_F_mat4)));
           
           subplot(2,2,1); hold on; imagesc(loc_F_mat, [low, high]); colorbar; title('F near saddle in YZ plane, X large');
           subplot(2,2,2); hold on; imagesc(loc_F_mat2, [low, high]); colorbar; title('F near saddle in YZ plane, X Middle-large');
           subplot(2,2,3); hold on; imagesc(loc_F_mat3, [low, high]); colorbar; title('F near saddle in YZ plane, X Middle-small');
           subplot(2,2,4); hold on; imagesc(loc_F_mat4, [low, high]); colorbar; title('F near saddle in YZ plane, X small');

           
           % Compute J AT the saddle point
           [J_AT_SADDLE , y_vec_AT_SADDLE]  = Jfun(x_vec, sigma, alpha, f_star, 0, soft_constrain_flag);
           
           %%return;

         
           
           
            % Start iterating and solving backwards. We know there IS a
            % critical point nearby. Why can't the numeric find it
            % ??????????
            f_star_is = f_star
            left_frac = f_vec(f_star_start)
            right_frac = f_vec(f_star_start+1)
            delta_v = zeros(1, 4); % The delta vector
            delta_v_approx = zeros(1,4); % The delta vector approximation
            a = zeros(1, 4); % A vector storing several integrals that we need
            
            % Compute the a's once and for all ! 
            a(2) = 2 * I * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, 0, 2);
            a(1) = 2 * I * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, 0, 2) + a(2);

            a(4) = 2 * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, 0, 3);
            a(3) = 2 * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, 0, 3) + a(4);
            a_is = a
            
            
            sigma_expansion = ( a(3).*a(2).^2 + a(4).*a(1).^2-2.*a(1).*a(2).*a(4)) ./ (alpha.^2 .* a(1).^2 .* 1 .* num_vars)
            sigma_expansion = sqrt(sigma_expansion)
%%%            sigma_expansion = sqrt(a(3).*a(2).^2 + a(4).*a(1)^2-2.*a(1).*a(2).*a(3)) ./ (alpha .* a(2) .* sqrt(num_vars))
            
            for f_star_ind = f_star_start:-1:1 %%%f_star_start-10 %%%%1
%                f_star_ind
                curr_frac = f_vec(f_star_ind)
                
                delta_v(4) = curr_frac-f_star; % update delta_f
                
                
                %%%%% Keep y and z zero and solve for the new x : 
%%%%                new_x = fzero('FUNC_for_x_only', x_vec(1), [],  f_vec(f_star_ind) , alpha, C_alpha, sigma)
%%%%                new_x_vec = x_vec; new_x_vec(1) = new_x;
                 
                x_vec
                sigma
                alpha
                f_vec(f_star_ind)
                0
                soft_constrain_flag
                new_x_vec = newtonSys('Jfun',x_vec, 0.000001, 0.000001, 250, [], sigma, alpha, f_vec(f_star_ind), 0, soft_constrain_flag);
                
                [FFF JJJ] = FInExponent( new_x_vec(1), new_x_vec(2), new_x_vec(3), sigma, alpha, C_alpha, f_vec(f_star_ind), 3, soft_constrain_flag);
                %should_be_positive_F = FFF 
                saddle_kept_frac_dist(f_star_ind) = sqrt(num_vars * det(J4) / (2*pi*det(JJJ))) * exp(-num_vars * FFF); 
                cur_prob = saddle_kept_frac_dist(f_star_ind)
               
                % Do work to debug the numerics using the Taylor expansion.
                % Note : the x_vec is given in the form of t, s and we want
                % to get back to y, z 
                delta_v(1:3) = new_x_vec; delta_v(1) = delta_v(1) - x_alpha;
                delta_v(2) = I * log(1 + delta_v(2));
                delta_v(3) = I * log(1 + delta_v(3)) - delta_v(2);
                delta_v_is = delta_v;
                
             %% quadl('SaddleHelperGaussian', C_alpha, 99, TOL, [], x_vec(1), sigma, x_vec(2), 2);                     
               
             
                % New - Break the integrals into several parts to attain
                % high accuracy !! 
%                 a(2) = 2 * I * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, 0, 2); 
%                 a(1) = 2 * I * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, 0, 2) + a(2);
%                            
%                 a(4) = 2 * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, 0, 3);
%                 a(3) = 2 * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, 0, 3) + a(4);    
%                 a_is = a;
               


                
                delta_v_approx(1) = a(1) * a(4) - a(2) * a(3);
                delta_v_approx(2) = a(1) * a(2);
                delta_v_approx(3) = -a(1)^2;
                delta_v_approx = (I*alpha ./ ( a(2)^2 * a(3) + a(1)^2 * a(4) - 2 * a(1) * a(2) * a(4))) .*  delta_v_approx; 
                delta_v_approx(4) = 1;
                delta_v_approx = delta_v(4) * delta_v_approx; 
                
                F_approx = 0.5 * (alpha * a(1) * delta_v(4))^2 / ( a(2)^2 * a(3) + a(1)^2 * a(4) - 2 * a(1) * a(2) * a(4))
                FFF_Numeric = FFF;
                
                % avoid extra work on zero probabilities
                if(cur_prob < TOL)
                    break;
                end
                
                x_vec = 2 .* new_x_vec - x_vec; % Change to the new starting point
                
                
            end
            
            
            x_vec = [x_alpha, 0, 0]';    % Init back the x_vec
            
            
            % Start iterating and solving forward
            for f_star_ind = f_star_start+1:length(f_vec) % f_star_start + 10 %%%%length(f_vec)
%                f_star_ind
                curr_frac = f_vec(f_star_ind)
                new_x_vec = newtonSys('Jfun',x_vec, 0.000001, 0.000001, 250, [], sigma, alpha, f_vec(f_star_ind), 0, soft_constrain_flag);
                [FFF JJJ] = FInExponent( new_x_vec(1), new_x_vec(2), new_x_vec(3), sigma, alpha, C_alpha, f_vec(f_star_ind), 3, soft_constrain_flag);
                FFF_Numeric = FFF; 
                saddle_kept_frac_dist(f_star_ind) = sqrt(num_vars * det(J4) / (2*pi*det(JJJ))) * exp(-num_vars * FFF); 
                cur_prob = saddle_kept_frac_dist(f_star_ind)
                
                 % avoid extra work on zero probabilities
                if(cur_prob < TOL)
                    break;
                end
                
                x_vec = 2 .* new_x_vec - x_vec; % Change to the new starting point
            end

        end
    end
else % Here we compare two versions of the measured correlations
    if(one_side_flag)
        if(rand_flag == UNIFORM)

        end
        if(rand_flag == LINEAR)

        end
        if(rand_flag == GAUSSIAN)

        end
    else  % here take two sides
        if(rand_flag == UNIFORM)

        end
        if(rand_flag == LINEAR)

        end
        if(rand_flag == GAUSSIAN)
            % Now solve the saddle point equations for two samples :
            % (complicated ....) 
            % Note : For now we do only two-sides with gaussian dists.
            % Initilize the X solution
            
            start_doing_two_samples = 999999
            x_vec = [x_alpha, x_alpha, 0, 0, 0]' %(x1, x2, y, z, w)
            
            
            %prev_x_vec; % A linear extrapulation trick to speed calculations
            % Here there are already problems (eventhough we give zeros as
            % input !!!!) 
            [THIS_F_SHOULD_BE_ZERO J6] = FInExponentTwoSampled( x_vec(1), x_vec(2), x_vec(3), x_vec(4), x_vec(5), sigma, alpha, C_alpha, f_star, 6,  soft_constrain_flag); 
            
            
            det_2nd_deriv_should_be_positive = det(J6);
            
            % find the index of the correct f_star
            f_star_start = sum(f_vec < f_star);

            if(f_star_start == 0)
                printf('Problem ! no f_vec left to f_star - Choose  a finer grid !!!\n')
            end
           
           % Compute J AT the saddle point - why is this needed ? 
           %%%%%[J_AT_SADDLE , y_vec_AT_SADDLE]  = Jfun(x_vec, sigma, alpha, f_star, 0, soft_constrain_flag)
           
            % Start iterating and solving backwards. We know there IS a
            % critical point nearby. Why can't the numeric find it
            % ??????????
            f_star_is = f_star
            left_frac = f_vec(f_star_start)
            right_frac = f_vec(f_star_start+1)
            delta_v = zeros(1, 6); % The delta vector
            delta_v_approx = zeros(1,6); % The delta vector approximation
            a = zeros(1, 4); % A vector storing several integrals that we need
            
            % Compute the a's once and for all ! why is this needed ?? 
            a(2) = 2 * I * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, 0, 2);
            a(1) = 2 * I * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, 0, 2) + a(2);

            a(4) = 2 * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_vec(1), sigma, 0, 3);
            a(3) = 2 * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_vec(1), sigma, 0, 3) + a(4);
            a_is = a
            
            % This is 'cheating' - we take sigma to be the true&sampled
            % ver. 
            sigma_expansion = ( a(3).*a(2).^2 + a(4).*a(1).^2-2.*a(1).*a(2).*a(4)) ./ (alpha.^2 .* a(1).^2 .* 1 .* num_vars)
            sigma_expansion = sqrt(sigma_expansion)
            
            for f_star_ind = f_star_start:-1:f_star_start-10 %%%%1
                curr_frac = f_vec(f_star_ind)
                
                delta_v(6) = curr_frac-f_star; % update delta_f
                
                
                %%%%% Keep y and z zero and solve for the new x : 
%%%%                new_x = fzero('FUNC_for_x_only', x_vec(1), [],  f_vec(f_star_ind) , alpha, C_alpha, sigma)
%%%%                new_x_vec = x_vec; new_x_vec(1) = new_x;
                 
%                 x_vec
%                 sigma
%                 alpha
%                 f_vec(f_star_ind)
%                 0
%                 soft_constrain_flag
%                dumb_new_x_vec = newtonSys('JfunTwoSampled',x_vec, 0.000001, 0.000001, 250, [], sigma, alpha, f_star, 0, soft_constrain_flag) % here we are supposed not to move !!!! 
                call_newton_sys = 99
                delta_f_is = f_vec(f_star_ind) - f_star
                new_x_vec = newtonSys('JfunTwoSampled',x_vec, 0.000001, 0.000001, 250, [], sigma, alpha, f_vec(f_star_ind), 0, soft_constrain_flag)
                finished_newton_sys = 999
                
                [FFF JJJ] = FInExponentTwoSampled( new_x_vec(1), new_x_vec(2), new_x_vec(3), new_x_vec(4), new_x_vec(5), sigma, alpha, C_alpha, f_vec(f_star_ind), 5, soft_constrain_flag);
                %should_be_positive_F = FFF
                FFF_Numeric = FFF;

                saddle_kept_frac_dist(f_star_ind) = sqrt(num_vars * det(J6) / (2*pi*det(JJJ))) * exp(-num_vars * FFF);
                cur_prob = saddle_kept_frac_dist(f_star_ind)

% % % % %                 % Do work to debug the numerics using the Taylor expansion.
% % % % %                 % Note : the x_vec is given in the form of t, s and we want
% % % % %                 % to get back to y, z 
% % % % %                 delta_v(1:3) = new_x_vec; delta_v(1) = delta_v(1) - x_alpha;
% % % % %                 delta_v(2) = I * log(1 + delta_v(2));
% % % % %                 delta_v(3) = I * log(1 + delta_v(3)) - delta_v(2);
% % % % %                 delta_v_is = delta_v;
% % % % % 
% % % % %                 
% % % % %                 delta_v_approx(1) = a(1) * a(4) - a(2) * a(3);
% % % % %                 delta_v_approx(2) = a(1) * a(2);
% % % % %                 delta_v_approx(3) = -a(1)^2;
% % % % %                 delta_v_approx = (I*alpha ./ ( a(2)^2 * a(3) + a(1)^2 * a(4) - 2 * a(1) * a(2) * a(4))) .*  delta_v_approx; 
% % % % %                 delta_v_approx(4) = 1;
% % % % %                 delta_v_approx = delta_v(4) * delta_v_approx; 
% % % % %                 
% % % % %                 F_approx = 0.5 * (alpha * a(1) * delta_v(4))^2 / ( a(2)^2 * a(3) + a(1)^2 * a(4) - 2 * a(1) * a(2) * a(4))
% % % % %                 FFF_Numeric = FFF
                
                % avoid extra work on zero probabilities
                if(cur_prob < TOL)
                    break;
                end
                
                x_vec = 2 .* new_x_vec - x_vec; % Change to the new starting point
                
                
            end
            
            
            x_vec = [x_alpha, x_alpha, 0, 0, 0]';    % Init back the x_vec
            
            
            % Start iterating and solving forward
            for f_star_ind = f_star_start+1:f_star_start + 10 %%%%length(f_vec)

                curr_frac = f_vec(f_star_ind)
                new_x_vec = newtonSys('JfunTwoSampled',x_vec, 0.000001, 0.000001, 250, [], sigma, alpha, f_vec(f_star_ind), 0, soft_constrain_flag);

                [FFF JJJ] = FInExponentTwoSampled( new_x_vec(1), new_x_vec(2), new_x_vec(3), new_x_vec(4), new_x_vec(5), sigma, alpha, C_alpha, f_vec(f_star_ind), 5, soft_constrain_flag);

                FFF_Numeric = FFF;

                saddle_kept_frac_dist(f_star_ind) = sqrt(num_vars * det(J6) / (2*pi*det(JJJ))) * exp(-num_vars * FFF);
                cur_prob = saddle_kept_frac_dist(f_star_ind)

                % avoid extra work on zero probabilities
                if(cur_prob < TOL)
                    break;
                end
                
                x_vec = 2 .* new_x_vec - x_vec; % Change to the new starting point
            end

            
            
        end
    end

end