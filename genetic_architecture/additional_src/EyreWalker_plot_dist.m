% Sample effect sizes and selection coefficient based on Eyre-Walker's
% distribution from Eyre-Walker, PNAS, 2010.
% http://www.pnas.org/content/107/suppl.1/1752.abstract
% The model is:
% s ~ Gamma(k, theta/k)
% beta = delta*s*(1+epsilon)
% epsilon ~ N(0, sigma_epsilon^2)
%
% Input:
% tau - power relating s and beta
% sigma_epsilon - noise in coupling of s and beta
% delta - constant relating s and beta (in the paper it is drawn uniformly as +/- 1 with prob. 1/2) 
% theta - scale parameter of the Gamma distribution
% k - shape parameter of the Gamma distribution
% N - effective population size (default is 10000) 
% n - number of loci to sample
% link_function - coupling between s and beta (default 'linear', also possible 'sigmoid')
% figs_dir - where to save output figure
%
% Output:
% beta - list of sampled effect sizes
% s - list of sampled selection coefficients
% f_vec - allele frequency vector
% V_vec - variance explained vector (as function of allele frequency) 
% f_var_bins - allele freqeuncy bins for histogram of variance explained
% f_var_hist - histogram of variance explained (as function of allele frequency)
%
function [beta s f_vec V_vec f_var_bins f_var_hist] = ...
    EyreWalker_plot_dist(tau, sigma_epsilon, delta, ... % parameters relating s and beta
    theta, k, N, n, link_function, figs_dir) % parametring specifying the gamma distribution for beta

if(~exist('N', 'var') || isempty(N))
    N = 10000; % effective population size of humans
end
if(~exist('link_function', 'var') || isempty(link_function)) % set default link function
    link_function = 'linear';
end
[beta s] = ...
    EyreWalker_dist_rnd(tau, sigma_epsilon, delta, ... % parameters relating s and beta
    theta, k, N, n, link_function); % call function to sample alleles 


mean_f_for_var_explained_given_s = mean_x_given_s(-s, N, 1); 



% Plot figure
title_str = ['Parameters: $\tau=' num2str(tau, 2) ', \sigma^2=' num2str(sigma_epsilon.^2,2) ...
    ', \delta=' num2str(delta,2) ', \theta=' num2str(theta, 4) ', k=' num2str(k, 2) '$'];
save_str = strdiff(title_str, '\'); save_str = strdiff(save_str, ',');
save_str = strdiff(save_str, '$'); save_str = strdiff(save_str, 'Parameters: '); save_str = strdiff(save_str, '^2');
save_str = strrep(save_str, ' ', '_'); save_str = strrep(save_str, '=', '_');


[beta_hist beta_bins] = hist_density(beta, 100, [], 0);
[s_hist s_bins] = hist_density(s, n/2500, [], 0);
[f_var_hist f_var_bins] = hist_density(mean_f_for_var_explained_given_s, n/2500, [], 0);

[f_var_hist_weighted f_var_bins_weighted] = weighted_hist(mean_f_for_var_explained_given_s, beta.^2, n/2500); 
f_var_hist_weighted = normalize_hist(f_var_bins_weighted, f_var_hist_weighted);

[f_var_hist_weighted_log f_var_bins_weighted_log] = weighted_hist(log(mean_f_for_var_explained_given_s), beta.^2, n/2500); 
f_var_hist_weighted_log = normalize_hist(f_var_bins_weighted_log, f_var_hist_weighted_log); % this is 'cheating' - normalize on the log scale. Just for display (does not yield a density function!)
f_var_bins_weighted_log = exp(f_var_bins_weighted_log); 


figure; subplot(2,2,4); % figure;
hold on; % xlabel('Effect size on trait (\beta)');
ylabel('Density');
bar(beta_bins, beta_hist, 'k'); % title(['Effect sizes distribution. ' title_str]);
x_lim = quantile(beta, 0.99); % don't display outliers
set(gca, 'Xlim', [0 x_lim]);
set(gca,'YDir','reverse'); % flip y axis
x_lim = get(gca, 'xlim');

subplot(2,2,1); % figure;
hold on; % ylabel('Selection Coefficient (s)');
xlabel('Density');
barh(s_bins, s_hist, 'k'); % title(['Selection coefficient distribution. ' title_str]);
y_lim = quantile(s, 0.99); % don't display outliers
set(gca, 'Ylim', [0 y_lim]);
set(gca,'XDir','reverse'); % flip x axis
y_lim = get(gca, 'ylim');

subplot(2,2,2); hold on;
% plot(beta, s, '.');
num_bins = 50;
x_bins = linspace(x_lim(1), x_lim(2), num_bins);
y_bins = linspace(y_lim(1), y_lim(2), num_bins);
hist2d_draw([s beta], y_bins, x_bins, [], [], [], 100 * n / num_bins^2, 0);
colormap ('gray'); colormap(flipud(colormap));
xlabel('Effect size on trait (\beta)'); ylabel('Selection coefficient (s)');
% set(gca, 'xlim', x_lim); set(gca, 'ylim', y_lim);

mtit(title_str, 'interpreter', 'latex', 'fontsize',14); % camroll(180)

res = 0.001; f_vec = res:res:0.5;
two_side_flag=1; % compute for minor allele frequency
V_vec=zeros(size(f_vec));
for i=1:length(s_bins) % compute variance explained distribution empirically
    cur_V = allele_freq_spectrum(f_vec, -s_bins(i), N, two_side_flag, [], 1);
    V_vec = V_vec + s_hist(i) .* cur_V;
end
V_vec = normalize_hist(f_vec, V_vec); % sum all variance to one


% V_vec_analytic = (theta./gamma(k)) .* (k./theta).^k .* (1+sigma_epsilon^2) .* gamma(2*tau+k) .* ...
%     (zeta_hurwitz(2*tau+k, f_vec + k ./ theta) - zeta_hurwitz(2*tau+k, 1 + k ./ theta)); % Formula using Hurwitz zeta function
% V_vec_analytic = normalize_hist(f_vec, vec2row(V_vec_analytic)); 

subplot(2,2,3); plot(f_vec, V_vec, 'linewidth', 2); xlim([0 0.5]); ylim([0 max(V_vec) .* 1.05]);
xlabel('Minor Allele Freq. (f)'); ylabel('Var. Explained (V)');
my_saveas(gcf, fullfile(figs_dir, ['EyreWalker_' link_function '_' save_str]), 'pdf'); % {'epsc', 'pdf', 'jpg'});


% New: Plot distribution of alleles according to s and distribution of
% variance explained according to s
s_var_hist = weighted_hist(s, beta.^2, [0 s_bins max(1, max(s) + 0.00001)]);
s_var_hist = normalize_hist(s_bins, s_var_hist(1:end-1));

mean_s = mean_hist(s_bins, s_hist);
mean_var_s = mean_hist(s_bins, s_var_hist);
median_s = median_hist(s_bins, s_hist);
median_var_s = median_hist(s_bins, s_var_hist);

for log_flag = 0:1
    for cumulative_flag = 0:1
        figure;
        switch log_flag
            case 0
                switch cumulative_flag
                    case 0 % plot density
                        log_str = '';
                        plot(s_bins, s_hist, 'linewidth', 2); hold on;
                        plot(s_bins, s_var_hist, 'r', 'linewidth', 2);
                    case 1 % plot cumulative
                        log_str = '_cum';
                        plot(s_bins, cumsum(s_hist) ./ sum(s_hist), 'linewidth', 2); hold on;
                        plot(s_bins, cumsum(s_var_hist) ./ sum(s_var_hist), 'r', 'linewidth', 2);
                end
            case 1
                switch cumulative_flag
                    case 0 % plot density
                        log_str = '_log_scale';
                        semilogx(s_bins, s_hist, 'linewidth', 2); hold on;
                        semilogx(s_bins, s_var_hist, 'r', 'linewidth', 2);
                    case 1
                        log_str = '_log_scale_cum';
                        semilogx(s_bins, cumsum(s_hist) ./ sum(s_hist), 'linewidth', 2); hold on;
                        semilogx(s_bins, cumsum(s_var_hist) ./ sum(s_var_hist), 'r', 'linewidth', 2);
                end
        end
        if(~cumulative_flag)
            ylim([0 max(s_var_hist) * 3]);
            cum_str = '';
            ylabel('Density');
        else
            cum_str = 'cum.';
            ylabel('Cumulative Prob.');
        end
        legend({'# alleles', 'var-explained'}, 4);
        xlabel('s');
        title(['Dist. of var. expl. and # alleles by selec. coeff. ' cum_str ' (mean s=' ...
            num2str(mean_s,3) ', mean var-s=' num2str(mean_var_s,3) ')']);
        my_saveas(gcf, fullfile(figs_dir, ['EyreWalker_' link_function ...
            '_var_explained_by_selective_coefficient_' save_str log_str]), 'pdf'); % {'epsc', 'pdf', 'jpg'});
    end % loop on cumulative flag
end % loop on log flag



for log_flag = 0:1 % New: Plot var. explained by allele frequency
    figure;
    switch log_flag
        case 0 % linear scale DAF
            log_str = ''; 
            plot(f_var_bins, f_var_hist); hold on;
            plot(f_var_bins_weighted, f_var_hist_weighted, 'r');
        case 1 % log-scale for DAF
            log_str = ' (log-scale)';
            semilogx(f_var_bins, f_var_hist); hold on;
            semilogx(f_var_bins_weighted_log, f_var_hist_weighted_log, 'r'); 
    end
    legend({'raw', 'weighted'});
    Title(['Distribution of Var. Explained by allele. freq.' log_str]); 
    xlabel(['Derived Allele Frequency' log_str]); ylabel('Var-Explained');
end

f_var_bins = f_var_bins_weighted_log; f_var_hist = f_var_hist_weighted_log; % choose output histogram 




