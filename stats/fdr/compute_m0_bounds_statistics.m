% Compute statistics for various m0 estimators: 
n=20; iters =10000;
mu = uniform_sum_reciprocal_mean(n, iters); mu=0.5.*mu(2:end);
mu_log = uniform_log_sum_reciprocal_mean(n, iters); mu_log = mu_log(2:end);
x_vec = 2:length(mu)+1;
figure; plot(x_vec, mu - (1./[2:length(mu)+1]')); title('mean minus 1/m_0');
figure; hold on; plot(x_vec, mu); plot(x_vec, mu_log, 'r'); 
title('mean '); legend('p_i', 'log');

figure; plot(log(x_vec),    log(mu - (1./[2:length(mu)+1]'))); title('log-log plot');

figure; plot(x_vec, mu.^2 ./ (mu-2)); title('corrected m_0');

figure; plot(x_vec, x_vec' - 1./mu); title('c vec');

