% Compute the peak of the delta distribution for the fraction when Ngenes goes to infinity.
% We also give the standard deviation according to the saddle-point, and the threshold x_alpha.
% Here dist_std is supposed to be the st.d. of the ORIGINAL Q, and nsamples
% should represent the st.d. of the additive noise
% Input:
% rand_flag - distributions of true correlations Q (GAUSSIAN, UNIFORM or LINEAR)
% one_side_flag - definition of 'top' genes (ONE_SIDE - only highest correlations, TWO_SIDES - highest ABSOLUTE correlations)
% true_corr_flag - type of overlap to compute: TRUE_AND_SAMPLED - overlap between true and noisy, TWO_SAMPLED - overlap between two noisy 
% dist_std - standard deviation of true distribution Q 
% nsamples - sample size
% alpha - fraction of top genes
% mu - mean of 
% prior - prior on distribution used to determine algorithm starting point 
%
% Output:
% inf_limit_frac - saddle-point approximation of mean overlap f
% inf_limit_std - saddle-point approximation of st.d. of overlap f
% x_alpha
%
function [inf_limit_frac, inf_limit_std, x_alpha] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, dist_std, nsamples, alpha, mu, prior)


global BOOTSTRAP NON_OVERLAP;
global UNIFORM GAUSSIAN LINEAR FROM_DATA mix_GAUSSIAN student_t;
global ONE_SIDE TWO_SIDES;
global TRUE_AND_SAMPLED TWO_SAMPLED;

TOL = 0.00000000001;

I = sqrt(-1);
inf_limit_std = 1; % Dummy for Matlab not to shout ! could be wrong !

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
    start_point=0;
    for i=1:length(prior)
        start_point=start_point+double(prior(i)*norminv(1-alpha, mu(i), dist_std(i)));
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
    start_point=0;
    for i=1:length(prior)
        start_point=start_point+double(prior(i)*norminv(1-0.5*alpha,mu(i),dist_std(i)));
    end    
end

if(rand_flag == mix_GAUSSIAN)
    end_x=1;
    another_end=1;
    while(another_end)
        [a, C_alpha, another_end]=myfsolve(mu, dist_std, prior, alpha, one_side_flag, end_x)
        end_x=end_x*1.5;
    end
    %     C_alpha = fsolve(@(c_a)find_c_alpha(mu,dist_std,prior,alpha,one_side_flag,c_a),start_point);
    sigma= 1 /sqrt(nsamples);
    noisy_std=sqrt(dist_std.^2+sigma.^2);
    another_end=1;
    while(another_end)
        [a, x_alpha, another_end] = myfsolve(mu, noisy_std, prior, alpha, one_side_flag, end_x)
        end_x=end_x*1.5;
    end
else
    sigma = 1 / (dist_std*sqrt(nsamples)); % st.d. is proportional to 1/sqrt(N)
end

if(true_corr_flag==TRUE_AND_SAMPLED)   % here we compare the true and measured correlations
    if(one_side_flag)
        if(rand_flag == UNIFORM)
            x_alpha = fzero('JointDensFracUniformDist', C_alpha + sigma*norminv(1-alpha),[],  sigma, alpha);
            inf_limit_frac = 1.0 - (1.0/(2*sqrt(3)*alpha))*(sigma*quad('normcdf', (x_alpha-sqrt(3))/sigma, (x_alpha-C_alpha)/sigma) );
        end
        if(rand_flag == LINEAR)
            x_alpha = fzero('JointDensFracLinearDist', C_alpha + sigma*norminv(1-alpha),[],  sigma, alpha);
            inf_limit_frac = 1.0 - (1.0/(alpha))*quad('lin_normcdf', C_alpha, sqrt(6), TOL, [], x_alpha, sigma, one_side_flag) ;
        end
        if(rand_flag == GAUSSIAN)
            x_alpha = C_alpha * sqrt(1+sigma.^2);
            inf_limit_frac = (1.0/alpha)*quadl('JointDensFrac', C_alpha, 100*(1+sqrt(sigma)),TOL, [], sigma, C_alpha, one_side_flag);
        end
        if(rand_flag == mix_GAUSSIAN)
            inf_limit_frac = (1.0/alpha)*quadl('MixJointDensFrac',C_alpha, 100*sum(noisy_std),TOL, [], ...
                sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std);
        end
    else  % here take two sides
        if(rand_flag == UNIFORM)
            x_alpha = fzero('JointDensFracUniformDist', C_alpha + sigma*norminv(1-alpha/2),[],  sigma, alpha/2);
            inf_limit_frac = 1.0 - (1.0/(sqrt(3)*alpha))*(sigma*quad('normcdf', (x_alpha-sqrt(3))/sigma, (x_alpha-C_alpha)/sigma) ) + ...
                (1.0/(sqrt(3)*alpha))*(sigma*quad('normcdf', (-x_alpha-sqrt(3))/sigma, (-x_alpha-C_alpha)/sigma) );
        end
        if(rand_flag == LINEAR)
            x_alpha = fzero('JointDensFracLinearDist', C_alpha + sigma*norminv(1-alpha/2),[],  sigma, alpha/2);
            inf_limit_frac = 1.0 - (1.0/(1*alpha))*quad('lin_normcdf', C_alpha, sqrt(6), TOL, [], x_alpha, sigma, one_side_flag) ;
        end
        if(rand_flag == GAUSSIAN)
            x_alpha = C_alpha * sqrt(1+sigma.^2);
            inf_limit_frac = (1.0/alpha)*quadl('JointDensFrac', C_alpha, 100*(1+sqrt(sigma)), TOL, [],...
                sigma, C_alpha, one_side_flag);
            % New ! (28.6.05) - Compute also the (approx.) standard
            % deviation analytically
            %%%  a = zeros(4);
            %%%  a(2) = 2 * I * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_alpha, sigma, 0, 2);
            %%%  a(1) = 2 * I * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_alpha, sigma, 0, 2) + a(2);
            
            %%%  a(4) = 2 * SaddleHelperGaussianQuadl(C_alpha, 99, TOL, x_alpha, sigma, 0, 3);
            %%%  a(3) = 2 * SaddleHelperGaussianQuadl(0, C_alpha, TOL, x_alpha, sigma, 0, 3) + a(4);
            
            %%%  inf_limit_std = sqrt( a(4) *(a(1)-a(2))^2 + (a(3) - a(4)) * a(2)^2 ) / (alpha * a(1));
            
            
            inf_limit_std = compute_inf_limit_std(rand_flag, one_side_flag, true_corr_flag, sigma, alpha,C_alpha, x_alpha, inf_limit_frac);
            
            
        end
        if(rand_flag == mix_GAUSSIAN)
            
            a = C_alpha; b=100*sum(noisy_std);
            subs_num = 10;
            dum_a = 1; dum_b = 1 + b-a;
            subs = 0:(log(dum_b)/(subs_num-1)):log(dum_b);
            subs = exp(subs)+a-1;
            inf_limit_frac = 0;
            
            for i=1:subs_num-1
                % Take the positive part
                inf_limit_frac = inf_limit_frac + quadl('MixJointDensFrac', subs(i), subs(i+1), TOL, [], ...
                    sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std);
                % Take the negative part
                inf_limit_frac = inf_limit_frac + quadl('MixJointDensFrac', -subs(i+1), -subs(i), TOL, [], ...
                    sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std);
                
            end
            
            inf_limit_frac = inf_limit_frac/ alpha;
            
            
            %             inf_limit_frac = (1.0/alpha)*(quadl('MixJointDensFrac', C_alpha, 100*sum(noisy_std), TOL, [], ...
            %                 sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std)+...
            %                 quadl('MixJointDensFrac', - 100*sum(noisy_std),-C_alpha, TOL, [], ...
            %                 sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std));
        end
    end
else % Here we compare two versions of the measured correlations
    if(one_side_flag)
        if(rand_flag == UNIFORM)
            x_alpha = fzero('JointDensFracUniformDist', C_alpha + sigma*norminv(1-alpha),[],  sigma, alpha);
            inf_limit_frac = (1.0/(alpha))*quad('uni_twosampled_normcdf', -sqrt(3), sqrt(3), TOL, [], x_alpha, sigma, one_side_flag) ;
        end
        if(rand_flag == LINEAR)
            x_alpha = fzero('JointDensFracLinearDist', C_alpha + sigma*norminv(1-alpha),[],  sigma, alpha);
            inf_limit_frac = (1.0/(alpha))*quad('lin_twosampled_normcdf', -sqrt(6), sqrt(6), TOL, [], x_alpha, sigma, one_side_flag) ;
        end
        if(rand_flag == GAUSSIAN)
            x_alpha = C_alpha * sqrt(1+sigma.^2);
            inf_limit_frac = (1.0/alpha)*quadl('gauss_twosampled_normcdf',-10*(1+sqrt(sigma)), 10*(1+sqrt(sigma)), TOL, [], C_alpha .* sqrt(1+sigma.^2), sigma, one_side_flag) ;
        end
        if(rand_flag == mix_GAUSSIAN)
            
            
            
            inf_limit_frac = (1.0/alpha)*(quadl('MixJointDensFrac', -100*sum(noisy_std), 100*sum(noisy_std), TOL, [],...
                sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std));
            
        end
    else  % here take two sides
        if(rand_flag == UNIFORM)
            x_alpha = fzero('JointDensFracUniformDist', C_alpha + sigma*norminv(1-alpha/2),[],  sigma, alpha/2);
            inf_limit_frac = (1.0/(1*alpha))*quad('uni_twosampled_normcdf', -sqrt(3), sqrt(3), TOL, [], x_alpha, sigma, one_side_flag) ;
        end
        if(rand_flag == LINEAR)
            x_alpha = fzero('JointDensFracLinearDist', C_alpha + sigma*norminv(1-alpha/2),[],  sigma, alpha/2);
            inf_limit_frac = (1.0/(1*alpha))*quad('lin_twosampled_normcdf', -sqrt(6), sqrt(6), TOL, [], x_alpha, sigma, one_side_flag) ;
        end
        if(rand_flag == GAUSSIAN)
            %             c_vec = [-999*sigma:sigma:999*sigma]; limit = 999*sigma
            %             brrr = gauss_twosampled_normcdf(c_vec, C_alpha, sigma, one_side_flag);
            %             figure; plot(c_vec, brrr);
            x_alpha = C_alpha * sqrt(1+sigma.^2);
            inf_limit_frac = (1.0/alpha)*quadl('gauss_twosampled_normcdf', -10*(1+sqrt(sigma)), -2*(1+sqrt(sigma)), TOL, [], C_alpha .* sqrt(1+sigma.^2), sigma, one_side_flag) +  ...
                (1.0/alpha)*quadl('gauss_twosampled_normcdf', -2*(1+sqrt(sigma)), 2*(1+sqrt(sigma)), TOL, [], C_alpha .* sqrt(1+sigma.^2), sigma, one_side_flag) + ...
                (1.0/alpha)*quadl('gauss_twosampled_normcdf', 2*(1+sqrt(sigma)), 10*(1+sqrt(sigma)), TOL, [], C_alpha .* sqrt(1+sigma.^2), sigma, one_side_flag);
            
            inf_limit_std = compute_inf_limit_std(rand_flag, one_side_flag, true_corr_flag, sigma, alpha,C_alpha, x_alpha, inf_limit_frac)
            
        end
        if(rand_flag == mix_GAUSSIAN)
            %             c_vec = [-999*sigma:sigma:999*sigma]; limit = 999*sigma
            %             brrr = gauss_twosampled_normcdf(c_vec, C_alpha, sigma, one_side_flag);
            %             figure; plot(c_vec, brrr);
            
            
            a = 0; b=100*sum(noisy_std);
            subs_num = 10;
            dum_a = 1; dum_b = 1 + b-a;
            subs = 0:(log(dum_b)/(subs_num-1)):log(dum_b);
            subs = exp(subs)+a-1;
            inf_limit_frac = 0;
            
            for i=1:subs_num-1
                % Take the positive part
                inf_limit_frac = inf_limit_frac + quadl('MixJointDensFrac', subs(i), subs(i+1), TOL, [], ...
                    sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std);
                % Take the negative part
                inf_limit_frac = inf_limit_frac + quadl('MixJointDensFrac', -subs(i+1), -subs(i), TOL, [], ...
                    sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std);
                
            end
            
            inf_limit_frac = inf_limit_frac/ alpha;
            
            % % %             inf_limit_frac = (1.0/alpha)*...
            % % %                 (quadl('MixJointDensFrac', - 100*sum(noisy_std), 100*sum(noisy_std), TOL, [], ...
            % % %                 sigma, x_alpha, one_side_flag,true_corr_flag,prior,mu,dist_std))
            
        end
    end
end
