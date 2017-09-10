%function test_regression
%regression_dir = '../../common_disease_model/docs/quasi_regression/';
regression_dir = '../../common_disease_model/docs/quasi_regression/tmp';

beta = [1 : 10]; % regression coefficients
n = 100000; % number of points
d = length(beta); % dimension
X = randn(10, 100000); % data
y =  X'*beta'; % outcome (dependent variable) 

beta_hat = (X*X')^(-1)*X*y; % regression 
beta_tilde = (X*y)./n; % quasi-regression 


x1_vec = X(1,:)' .* y; 
beta_tilde_sqr = sum(x1_vec)^2 - x1_vec'*x1_vec
beta_tilde_sqr_mean = beta_tilde_sqr / (n*(n-1))

    
% Take a non-linear model
beta12 = 1; 
YY = X(1,:) * beta(1) + X(2,:) * beta(2) + X(1,:).*X(2,:)*beta12;
beta_hat = (X(1:2,:)*X(1:2,:)')^(-1)*X(1:2,:)*YY' % regression 
beta_tilde = (X(1:2,:)*YY')./n % quasi-regression 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on n and d and estimate variance 
n_vec = 100:100:1000; num_n = length(n_vec); % vector of number of data points
d_vec = [500 5000]; num_d = length(d_vec); % vector of dimension 
iters = 50;

time_regression = cputime;
[mean_beta std_beta true_bias_beta estimated_bias_beta true_std_beta] = simulate_quasi_regression(n_vec, d_vec, iters);
time_regression = cputime - time_regression 

figure; imagesc(true_bias_beta); colorbar; title('Bias in \beta quasi-regression estimator');
xlabel('d'); ylabel('n');

figure; imagesc(std_beta); colorbar; title('Variance in \beta quasi-regression estimator');
xlabel('d'); ylabel('n');


legend_vec = [repmat('d=', num_d, 1) num2str(d_vec')];
figure; plot(n_vec, mean_beta, 'linewidth', 2); hold on; 
legend(legend_vec); 
xlabel('n'); ylabel('$\mu(\tilde{V_1})$', 'interpreter', 'latex'); 
title('Mean of $\tilde{V_1}$ as function of sample size. (un-adjusted)', 'interpreter', 'latex'); 
my_saveas(gcf, fullfile(regression_dir, 'figs', 'bias_linear_uncorrected'), {'epsc', 'pdf', 'fig'});

figure; plot(n_vec, mean_beta-estimated_bias_beta, 'linewidth', 2);
legend(legend_vec); 
xlabel('n'); ylabel('$\mu(\tilde{V_{1,BC}})$', 'interpreter', 'latex'); title('Mean of $\tilde{V}_{1,BC}$ as function of sample size. bias-corrected', 'interpreter', 'latex'); 
my_saveas(gcf, fullfile(regression_dir, 'figs', 'bias_linear'), {'epsc', 'pdf', 'fig'});

figure; plot(n_vec, std_beta, 'linewidth', 2); hold on; 
legend(legend_vec); 
plot(n_vec, true_std_V1, '--', 'linewidth', 2); 
xlabel('n'); ylabel('$\sigma(\tilde{V_1})$', 'interpreter', 'latex'); 
title('Standard Deviation of $\tilde{V_1}$ as function of sample size. (un-adjusted)', 'interpreter', 'latex'); 
my_saveas(gcf, fullfile(regression_dir, 'figs', 'variance_linear'), {'epsc', 'pdf', 'fig'});



% figure; plot(n_vec, mean_beta(:, 1), '*'); hold on; plot(n_vec, mean_beta(:, 1), 'b');
% plot(n_vec, mean_beta(:, 2), '*r'); hold on; plot(n_vec, mean_beta(:, 2), 'r');
% plot(n_vec, mean_beta(:, 1)-estimated_bias_beta(:,1), '*g'); hold on; plot(n_vec, mean_beta(:, 1)-estimated_bias_beta(:,1), 'g', 'linewidth', 3);
% plot(n_vec, mean_beta(:, 2)-estimated_bias_beta(:,2), '*m'); hold on; plot(n_vec, mean_beta(:, 2)-estimated_bias_beta(:,2), 'm', 'linewidth', 3);



figure; plot(n_vec, std_beta(:, 1), '*'); hold on; plot(n_vec, std_beta(:, 1), 'b');
plot(n_vec, std_beta(:, 2), 'r*'); hold on; plot(n_vec, std_beta(:, 2), 'r');
xlabel('n'); ylabel('$\sigma(\tilde{V_1})$', 'interpreter', 'latex'); title('Error in estimating linear variance as function of sample size'); 










