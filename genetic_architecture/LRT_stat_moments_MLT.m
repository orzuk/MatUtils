% Compute the mean and variance of the LoglikelihoodRaTio statistic between the MLT and LT model
%
% Input:
% N - total # of liabilities
% K - # liabilities needed to exceed threshold
% mu - prevalence
% h_x - heritability of x_1 assuming LT model
% true_model_str - model generating the data (MLT or LT) 
%
% Output:
%
% LRT_mu - mean of the log-likelihood-ratio test statistic
% LRT_var - variance of the log-likelihood-ratio test statistic
%
function [LRT_mu LRT_var] = LRT_stat_moments_MLT(N, K, mu, h_x, true_model_str)

epsilon=10^(-10); %AssignGeneralConstants;

if(~exist('true_model_str', 'var') || isempty(true_model_str))
    true_model_str = 'MLT';
end

x_mu = norminv(1-mu);
if(length(h_x) == 1) % we allow to already enter the heritability on the MLT model (one liability) 
    h_x(2) = heritability_scale_change_MLT(N*h_x, K, N, mu, 'MLT'); % Get the heritabiity of each liability in MLT model
end
mu(2) = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu)), 0, 1); % find mu_l that keeps the prevalence
x_mu(2) = norminv(1-mu(2)); % set threshold
N = [1 N]; 
switch true_model_str
    case 'MLT'
        true_ind = 2;
    case 'LT'
        true_ind = 1;
end
x_min = -5; x_max = 5; % integral limits

% % n=200^2;
% % [x1 x2] = meshgrid(norminv(((1:sqrt(n))-0.5) ./ sqrt(n)));
% % LRT_mu_mat =  ... % matrix of the LRT along the grid 
% %     z_expected_given_two_x_MLT(x1, x2, N(true_ind), h_x(true_ind), mu(true_ind), x_mu(true_ind)) .* ...
% %     ( log(max(epsilon^2,z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)))) - ...
% %     log(max(epsilon^2,z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)))) ) + ...
% %     (1-z_expected_given_two_x_MLT(x1, x2, N(true_ind), h_x(true_ind), mu(true_ind), x_mu(true_ind))) .* ...
% %     ( log(max(epsilon^2,1-z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)))) - ...
% %     log(max(epsilon^2,1-z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)))) );
% % figure; imagesc(LRT_mu_mat); colorbar;
% % figure; surf(x1, x2, LRT_mu_mat);
% % mean_LRT = mean(LRT_mu_mat(:))
% % var_LRT = var(LRT_mu_mat(:))

% 1 is the true model
LRT_mu = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
    ( z_expected_given_two_x_MLT(x1, x2, N(true_ind), h_x(true_ind), mu(true_ind), x_mu(true_ind)) .* ...
    ( log(max(epsilon^2,z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)))) - ...
    log(max(epsilon^2,z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)))) ) + ...
    (1-z_expected_given_two_x_MLT(x1, x2, N(true_ind), h_x(true_ind), mu(true_ind), x_mu(true_ind))) .* ...
    ( log(max(epsilon^2,1-z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)))) - ...
    log(max(epsilon^2,1-z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)))) ) ), ...
    x_min, x_max, x_min, x_max);

LRT_var = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
    ( z_expected_given_two_x_MLT(x1, x2, N(true_ind), h_x(true_ind), mu(true_ind), x_mu(true_ind)) .* ...
    ( log(max(epsilon^2,z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)))) - ...
    log(max(epsilon^2,z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)))) ).^2 + ...
    (1-z_expected_given_two_x_MLT(x1, x2, N(true_ind), h_x(true_ind), mu(true_ind), x_mu(true_ind))) .* ...
    ( log(max(epsilon^2,1-z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)))) - ...
    log(max(epsilon^2,1-z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)))) ).^2 ), ...
    x_min, x_max, x_min, x_max);
LRT_var = LRT_var - LRT_mu^2; 


% z_expected_LT = z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu);
% z_expected_MLT = z_expected_given_two_x_MLT(x1, x2, N, h_x_MLT, mu_l, x_mu_l);
%
% LRT_mu = z_expected_LT .* (log(z_expected_MLT)-log(z_expected_LT)) + ...
%     (1-z_expected_LT) * (log(1-z_expected_MLT)-log(1-z_expected_LT));
% LRT_var = LRT_mu*(1-LRT_mu); % This is wrong !!!!
