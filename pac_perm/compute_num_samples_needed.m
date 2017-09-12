% A lot of calculations, to compute the numbers of samples one needs in order to get
% a certain overlap 1-\eps with confidence 1-\delta
% We cannot solve the desired equarions analytically, therefore the
% bisection method is used for the optimization
function [num_samples_needed std_needed] = compute_num_samples_needed(rand_flag, one_side_flag, true_corr_flag, sigma, alpha, epsilon, delta, Ngenes, miu, prior)
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
TOL = 0.000001;
I = sqrt(-1);

TWO =2;

% Solve using the bisection method
sigma_left = 0.000001;
sigma_right = 10;

epsilon
delta

while(sigma_right-sigma_left > TOL)

    % Get the new sigma
    sigma_mid = (sigma_left + sigma_right)/2;
    % Compute f_star
    [f_star f_sigma x_alpha ] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, sigma, 1/(sigma_mid^2), alpha,miu,prior)  % We use here the approximation of Fisher
    
%    Ngenes
    f_sigma = f_sigma / sqrt(Ngenes);
    
    % Compute now the probability
    
%     is_legal = (1-epsilon-f_star)/f_sigma
%     epsilon
%     f_star 
%     f_sigma
     prob_f_above_thresh = 1 - normcdf((1-epsilon-f_star)/f_sigma);
    
    if(prob_f_above_thresh > 1 - delta)
        sigma_left = sigma_mid;
    else
        sigma_right = sigma_mid;
    end
        
    sigma_right
    sigma_left
    diff_sigma = sigma_right - sigma_left
    
end
    

% Return the new sigma
std_needed = sigma_mid
num_samples_needed = 1/sigma_mid^2 + 3  % Based on Fisher ..
