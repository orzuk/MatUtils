% Example of simulating a population with genotypes, IBD-blocks and phenotypes
% and estimating narrow-sense heritability for this population
plot_flag=0; num_bins = 100; % for plotting

%%%%%%%%%%%%% Set population demographic parameters (similar but higher IBD compared to paper's example) %%%%%%%
num_founders = 10; IBD_mean = 1/num_founders; % determines average IBD sharing
num_generations = 10;
num_individuals = 1000;
f_vec = mat2vec(repmat([0.15 0.31 0.21 0.21 0.5 0.17 0.2 0.13 0.31 0.41 0.2 0.42 0.2 ...
    0.17 0.15 0.15 0.4 0.5 0.2 0.2 0.22 0.2 0.17 0.2 0.11], 2, 1)); % Take a range of different allele frequencies

%%%%%%%%%%%%% Set LP model parameters from the paper's example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trait_type = 'quantitative';
num_pathways = 4; % Number of pathways in LP model
h_pathway = 0.5;
h_shared_env = 0.5 * (1-h_pathway);
[~, ~, ~, ~, ~, h_all_true] = ...
    compute_k_of_N_gaussian_statistics(0, 1, h_pathway, h_shared_env, [], ...
    'MAX', [], num_pathways, 10, 'numeric', 0, {'ACE'}); % Compute h_all (narrow sense heritability) for quantitative trait
num_snps = length(f_vec); num_causal_snps = num_snps/2;
beta_vec = ones(num_snps,1); beta_vec(1:num_snps/2) = 0; % only first half of genotypes contribute to phenotype 

%%%%%%%%%%%%% Simulate SNPs and estimate heritability. Repeat many iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iters=100; h_all_estimated = zeros(iters, 1); 
true_IBD_mat=[]; input_founder_identities = []; input_poiss_points = []; input_SNP_founder_identity_mat = [];
for i=1:iters
    if(mod(i,10)==0)
        run_iter = i
    end
    [true_IBD_mat, ~, SNP_mat, ...  % Simulate IBD matrix and many SNP matrices
        founder_SNP_mat, input_SNP_founder_identity_mat, input_founder_identities, input_poiss_points] = ...
        simulate_IBD_blocks_sharing(num_founders, num_generations, num_individuals, f_vec, ...
        true_IBD_mat, input_founder_identities, input_poiss_points, input_SNP_founder_identity_mat);
    
    LP_phenotype_vec = simulate_LP_phenotypes(SNP_mat, f_vec, num_pathways, h_pathway, h_shared_env, trait_type); % Simulate phenotype
    [h_all_estimated(i) phenotype_corr_by_IBD_vec] = estimate_heritability_by_local_regression(true_IBD_mat, LP_phenotype_vec, num_bins, plot_flag); % estimate heritability

    update_all_iters_phenotype_corr_by_IBD_vec; % update counts for all iterations
end

%%%%%%%% Plot results: mean and st.d. of estimated heritability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; boxplot(h_all_estimated); 
hold on; plot(1, h_all_true, 'k*');
ylabel('Estimated narrow-sense heritability');
title(['Estimated h^2 by local regression in ' num2str(iters) ' simulations']);
legend('true h_{all}^2');


%%%%%%%% Plot average regression line (similar to Figure 2 in pnas paper) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First average over all iterations
all_iters_phenotype_corr_by_IBD_vec.corr_mean = all_iters_phenotype_corr_by_IBD_vec.corr_mean ./ ...
    all_iters_phenotype_corr_by_IBD_vec.corr_counts;
all_iters_phenotype_corr_by_IBD_vec.corr_std = ...
sqrt( (all_iters_phenotype_corr_by_IBD_vec.corr_std - all_iters_phenotype_corr_by_IBD_vec.corr_mean.^2) ./ ...
all_iters_phenotype_corr_by_IBD_vec.corr_counts ) ./ sqrt(all_iters_phenotype_corr_by_IBD_vec.corr_counts);

figure; hold on;
x_vec = linspace(0, 1, 100); % why estimated_IBD_mat goes beyond one???
plot_beta = median(h_all_estimated) / (1-IBD_mean)^2; % Compute median beta
plot_y0 = -plot_beta .* IBD_mean; % make sure regression line passes at (k0, 0)

% Normalize slope and multiply by constant to get an estimator for heritability
plot_mean_h_slope = plot_beta*(1-IBD_mean)^2;
[~, plot_beta_std] = median_mad(h_all_estimated); plot_beta_std = plot_beta_std / (1-IBD_mean);  % compute MAD robust std for beta 

populated_bins = find(all_iters_phenotype_corr_by_IBD_vec.corr_counts >= 1000);
errorbar(all_iters_phenotype_corr_by_IBD_vec.IBD_bins(populated_bins), all_iters_phenotype_corr_by_IBD_vec.corr_mean(populated_bins), ...
    all_iters_phenotype_corr_by_IBD_vec.corr_std(populated_bins), ...
    'b', 'linewidth', 2); % plot averages of phenotypic correlations vs. IBD sharing 
plot(x_vec, plot_beta .* x_vec + plot_y0, 'r', 'linewidth', 2); % plot mean regression line

xlabel('IBD sharing'); ylabel('\rho(IBD)');
legend('\rho empirical', 'slope at k_0', 2); legend boxoff;
y_lim = get(gca, 'ylim');
line([IBD_mean IBD_mean], [y_lim(1)  0], ...
    'color', 'k', 'linestyle', '--');
line([0 IBD_mean], [0 0], ...
    'color', 'k', 'linestyle', '--');
text(IBD_mean, y_lim(1)*1.1, 'k_0', 'fontsize', 14);
xlim([0,1]);
slope_str = ['$h_{all}^2=' num2str(h_all_true,3) '$;  $h_{slope(mean-IBD)}^2 =' num2str(plot_mean_h_slope,3) ...
    ' \pm ' num2str(plot_beta_std,3) '$'];
text(IBD_mean*2, IBD_mean*0, slope_str, 'fontweight', 'bold', 'Interpreter','latex');


