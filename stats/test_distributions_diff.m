% A script for testing the distributions_diff function
% x - data sampled from a mixture-of-gaussians model
% x1 - the first gaussian (known - assumed to be neutral distribution)
% x2 - the second gaussian (unknown - assumed to be the 'selected' distribution)
LEFT = 0; RIGHT = 1;
n = 250000; mu = [0 2]; sigma = [1 1.7]; P = [0.8 0.2];   fdr_q = 0.05; eps = 0.5; num_bins = 500;

% x2 = randn(n, 1) .* sigma(1) + mu(1);
% x1 = randn(n, 1) .* sigma(2) + mu(2);
% x1 = MixtureOfGaussiansSimulateData(P, mu, sigma, n)'; % x2 = randn(100000,1) .* 0.5 - 2;
[x G] = MixtureOfGaussiansSimulateData(P, mu, sigma, n); x = x';
x1 = x(G == 1); x2 = x(G == 2);

[h1 bins1] = hist(x1, num_bins); [h2 bins2] = hist(x2, num_bins);
Threshold = TwoGaussiansThreshold( mu(1), sigma(1), mu(2), sigma(2))
IntegralsDiff = TwoGaussiansDiff( Threshold,  mu(1), sigma(1), mu(2), sigma(2), 1) % take right tail integral

figure; hold on; [p b] = hist_density(x, num_bins); % determine bins
[p1 b1] = hist_density(x1, b, 'g', 1, P(1)); [p2 b2] = hist_density(x2, b, 'r', 1, P(2));
legend('F_G', 'F_0', 'F_1'); xlabel('val.'); ylabel('freq.');
title('Comparison of separation methods'); 
my_saveas(gcf, 'two_gaussians.eps', 'epsc'); 

% Try Gaussian smoothing
gauss_smooth_width = 2; p_g = gauss_smooth(p, gauss_smooth_width); w= 10; 
b_g = [(b(1:10) + b(1)-b(11)) b (b(end-9:end) + b(end) - b(end-10))];
figure; hold on; plot(b,p, 'b'); plot(b_g,p_g, 'r'); title('Gaussian Smoothing Test');
legend(['orig. data (\mu = ' num2str(mean_hist(b,p)) ')'], ['smoothed data (\mu = ' num2str(mean_hist(b_g,p_g)) ')']);  


% True, FDR, conservative, liberal, MoG, epsilon
seperation_method_vec = {'FDR', 'conservative', 'liberal', 'MoG', ['\epsilon = ' num2str(eps)]};
figure; hold on; 
subplot(2,3,1); hold on; % plot true mixture 
% [p b] = hist_density(x, num_bins); % determine bins
% right_dist = hist_density(x(G == 2), num_bins,[],0); % compute the right Gaussian (don't plot)
% right_diff = IntegralsDiff; left_diff = IntegralsDiff;
% plot(b1, p1 .* (1-right_diff(1)), 'g'); % plot the fitted curve        [p1 b1] = hist_density(x1, b, 'g', 1, P(1));
% plot(b, right_dist .* right_diff(1), 'r'); % plot the fitted curve    [p2 b2] = hist_density(right_dist{i}, b, 'r', 1, P(2));
% legend('F_G', 'F_0', 'F_1');
% title(['True: Right diff = ' num2str(right_diff(1)) ', Left diff = ' num2str(left_diff(1)) ]);

fig_handle = MixtureOfGaussiansDraw1dGaussians(x, P, mu, sigma, ...
    {'', ''}, {'F_0', 'F_1', 'F_G'}, 'grb', [], num_bins, 0,0);
title(['True: Right diff = ' num2str(P(2)) ', Left diff = ' num2str(P(1)) ]);
right_dist = {};
diff_param = [fdr_q, fdr_q, fdr_q, fdr_q, eps];
for i=1:5
    [right_diff(i) right_cut(i) right_x{i} right_dist{i}] = distributions_diff(b, p, b1, p1, i, RIGHT, diff_param(i), 1, 1);
    [left_diff(i) left_cut(i) left_x{i} left_dist{i}] = distributions_diff(b, p, b1, p1, i, LEFT, diff_param(i), 1, 1);
    
    subplot(2,3,i+1); hold on; [p b] = hist_density(x, num_bins); % determine bins
    plot(b1, p1 .* (1-right_diff(i)), 'g'); % plot the fitted curve        [p1 b1] = hist_density(x1, b, 'g', 1, P(1));
    plot(b, right_dist{i} .* right_diff(i), 'r'); % plot the fitted curve    [p2 b2] = hist_density(right_dist{i}, b, 'r', 1, P(2));
    legend('F_G', 'F_0', 'F_1');
    title([seperation_method_vec{i} '. Right diff = ' num2str(right_diff(i)) ', Left diff = ' num2str(left_diff(i)) ]);
end
right_diff_is = right_diff
my_saveas(gcf, 'two_gaussians', 'epsc'); 

%    ', Right FDR = ' num2str(right_FDR) ', Left FDR = ' num2str(left_FDR) ', (fdr q = ' num2str(fdr_q) ')']);

labels_vec = {'cons.', 'freq.'};
legends_vec = {'',  'AR', 'Selected', 'Genomic',};
MixtureOfGaussiansDraw1dGaussians(x, P, mu, sigma, ...
    labels_vec, legends_vec, [], [], num_bins)


% Test the new randomization function
% figure; hold on;
% x1 = MixtureOfGaussiansSimulateData(P, mu, sigma, n)'; % randn(n,1);
% [p1 y1] = hist_density(x1, 10*num_bins, [], 0);
% x2 = distrnd(y1, p1, n, 1);
% [p2 y2] = hist_density(x2, num_bins, 'r', 0);

