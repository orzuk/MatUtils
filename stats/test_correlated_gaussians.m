% Generate correlated gaussians. Then perform de-convulution
% function test_correlated_gaussians()
test_gaussian = 0;
test_LD_association = 1;

if(test_gaussian)
    n_causal_vars = 10; % number of different causal SNPs (Gaussian variables)
    n_vars = 1000; % number of different total SNPs variables
    n_samples = 1000; % sample size (determined variance)
    beta_vec = zeros(n_vars,1); beta_vec(1:n_causal_vars) = randn(n_causal_vars,1); % original effect sizes
    C = randn(n_vars); C = C'*C; isposdef(C) % generate a correlation matrix
    beta_sigma_vec = diag(C); % variance in original beta's
    mu_vec = C * beta_vec; % means of different Gaussians
    C = C ./ (n_samples); % variance and covariance are proportional to 1/n_samples
    
    iters = 2000; % number of realizations for each variable
    
    mu_hat = simulate_correlated_gaussians(iters, n_vars, C, mu_vec); % these are the observed effect sizes we see in GWAS
    
    mu_hat_mean = vec2column(mean(mu_hat)); % mean observed effect size (should be very close to true effect size)
    C_hat = cov(mu_hat);
    
    figure; plot(C(:), C_hat(:), '.'); xlabel('original C_{ij}'); ylabel('observed C_{ij}'); % see that correlations were estimated correctly
    figure; plot(mu_vec, mu_hat_mean, '.'); xlabel('original \mu (with LD)'); ylabel('observed \mu (with LD)');
    
    %C_inv = inv(C); % get inverse matrix
    
    beta_orig_hat = C\mu_hat'; beta_orig_hat_mean = mean(beta_orig_hat,2);
    figure; plot(beta_vec, beta_orig_hat_mean, '.');
    xlabel('original \beta'); ylabel('observed \beta');
    
    
    true_var_explained = sum(beta_vec.^2)
    estimated_var_explained = sum(mu_hat.^2, 2)
    
    figure; hold on; hist_density(estimated_var_explained, 100);
    plot(true_var_explained, 0, 'r*', 'linewidth', 4);
end % test gaussians

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(test_LD_association) % New: Test genotypes associations
    trait_mode = [];
    f_vec = [0.1 0.1 0.3 0.3 0.4 0.4 0.2 0.2 0.15 0.15];
%    f_vec(:) = 0.1; % simplify: all are equal 
    beta_vec = [0 0.4 0.5 0.7 0 0 0 0 0 0]; % effect size
    n_samples = 1000000;
    num_snps = length(f_vec);
    C = posdefrnd(num_snps); % C = eye(num_snps); % correlation matrix
    SNP_CORR_mat = corr_mat_to_IBD_mat(C, f_vec); % This is the r^2 matrix 
    
    SNP_beta_mat = (SNP_CORR_mat .* ( sqrt(f_vec .* (1-f_vec))' * (1./sqrt(f_vec .* (1-f_vec))) ))'; 
    figure; plot(SNP_CORR_mat(:), SNP_beta_mat(:), '.')
    
    
    figure; imagesc(C); colorbar; % get correlation matrix
    figure; imagesc(SNP_CORR_mat); colorbar; % get IBD-sharing matrix
    figure; plot(C(:), SNP_CORR_mat(:), '.'); % see the relation between IBD and correlation
    xlabel('Gaussian correlation'); ylabel('Bernoulli correlation');
    
    h_vec = beta_vec.^2 .* f_vec .* (1-f_vec);
    h = sum(h_vec); % determine total heritability of trait
%    h_cov = sum(sum((beta_vec' * beta_vec) .* SNP_CORR_mat .* ...
%        sqrt( (f_vec .* (1-f_vec))' * (f_vec .* (1-f_vec)) ) ))
    h_cov = sum(sum((beta_vec' * beta_vec) .* SNP_beta_mat .* ...
        sqrt( (f_vec .* (1-f_vec))' * (f_vec .* (1-f_vec)) ) ))
    
    
    
    % [SNP_mat] = ... % corr_mat
    %     simulate_IBD_genotypes(K, f_vec, iters); % simulate
    % SNP_mat = sampleDichGauss01(f_vec,C, n_samples);
    % function
    % SNP_mat = simulate_correlated_gaussians(n_samples, length(f_vec), C, zeros(length(f_vec),1));
    % SNP_mat = SNP_mat > repmat(x_f_vec, n_samples, 1);
    
    SNP_mat = simulate_correlated_bernoulli(n_samples, length(f_vec), C, f_vec); % simulate genotypes
    SNP_mat_f = SNP_mat - repmat(f_vec, n_samples, 1);
    phenotype_vec = sum(SNP_mat_f .* repmat(beta_vec, n_samples, 1),2) + 1*randn(n_samples, 1) .* sqrt(1-h_cov);
    
    %SNP_mat = simulate_correlated_bernoulli(n_samples, 2, C, f_vec);
    
    C_empirical = corr(SNP_mat);
    figure; imagesc(corr(SNP_mat)); colorbar;
    figure; hold on; plot(mat2vec(SNP_CORR_mat), C_empirical(:), '.');
    plot(0:1, 0:1, 'r');
    xlabel('True correlations'); ylabel('Empirical correlations');
    f_vec_empirical = mean(SNP_mat);
    figure; hold on; plot(f_vec, f_vec_empirical, '.');
    plot(0:1, 0:1, 'r');
    xlabel('True frequencies'); ylabel('Empirical frequencies');
    
    empirical_beta_vec = zeros(num_snps, 1);
    for i=1:num_snps% Estiamte beta's
        empirical_beta_vec(i) = corr(SNP_mat(:,i), phenotype_vec);
    end
    empirical_beta_vec = empirical_beta_vec ./ vec2column( sqrt(f_vec .* (1-f_vec)) );
    reconstructed_beta_vec = empirical_beta_vec' / SNP_beta_mat'; % sqrtm(SNP_CORR_mat);
    figure; hold on;     plot(beta_vec, empirical_beta_vec, 'o');
    plot(beta_vec, reconstructed_beta_vec, 'r*');
    plot(0:1, 0:1, 'r'); 
    xlabel('True \beta'); ylabel('Empirical \beta'); legend('empirical', 'de-convolved');     
    
    effective_beta_vec = SNP_beta_mat*vec2column(beta_vec);
    figure; hold on; plot(effective_beta_vec, empirical_beta_vec, '.');
    xlabel('effective \beta'); ylabel('observed \beta');
    plot(0:1, 0:1, 'r'); 
    
    
    
    
    % SNP_CORR_mat,n_samples);
    % [contigency_table, genotype_phenotype_prod, genotype_sqr] = ...
    %     simulate_genotype_phenotype(p_vec, n_samples, [], ...
    %     iters, 'quantitative', trait_mode, 1)
    % phenotype_vec = contigency_table(:,1); % get the average phenotypes
    %
    %
    
    % Perform deconvulution
    
end % test LD association
