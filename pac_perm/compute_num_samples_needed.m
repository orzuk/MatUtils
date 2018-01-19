% Compute the numbers of samples needed in order to get
% a certain overlap 1-\eps with confidence 1-\delta
% The bisection method is used for the optimization
%
% Input:
% rand_flag - distributions of true correlations Q (GAUSSIAN, UNIFORM or LINEAR)
% one_side_flag - definition of 'top' genes (ONE_SIDE - only highest correlations, TWO_SIDES - highest ABSOLUTE correlations)
% true_corr_flag - type of overlap to compute: TRUE_AND_SAMPLED - overlap between true and noisy, TWO_SAMPLED - overlap between two noisy 
% sigma - standard deviation of true distribution Q 
% alpha - fraction of top genes
% epsilon - error level
% delta - confidence level
% Ngenes - total # of genes
% mu - mean of 
% prior - prior on distribution used to determine algorithm starting point 
%
% Output:
% num_samples_needed - sample size required 
% std_needed - 
%
function [num_samples_needed, std_needed] = compute_num_samples_needed(rand_flag, one_side_flag, true_corr_flag, sigma, alpha, epsilon, delta, Ngenes, mu, prior)
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
TOL = 0.000001;
I = sqrt(-1);

TWO =2;

% Solve using the bisection method
sigma_left = 0.000001;
sigma_right = 10;

while(sigma_right-sigma_left > TOL)
    % Get the new sigma
    sigma_mid = (sigma_left + sigma_right)/2;
    % Compute f_star
    [f_star, f_sigma, x_alpha] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, sigma, 1/(sigma_mid^2), alpha,mu,prior)  % We use here the approximation of Fisher    
    f_sigma = f_sigma / sqrt(Ngenes);
    
    % Compute now the probability    
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
num_samples_needed = 1/sigma_mid^2 + 3  % Based on Fisher transformation
