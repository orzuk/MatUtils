% Estimator should be unbiased !!!

% function test_estimate_IBD_sharing_from_snps
AssignGeneralConstants;
snp_type = 'binary';
iters = 1;
max_generations = 1; % how large a family to simulate
IBD_mean = 0.02; % average relatedness
IBD_std = 0.1; % standard deviation of relatedness
do_one_family=0; % run one family with IBD=1/2, 1/4 etc. to see that it's esitmated correctly
do_slope_analysis = 0; % run many families/pairs of individuals, assume we know IBD and regress
do_all_pairs = 0; % slope analysis for all pairs of individuals
analyze_REML_results = 0; % see what the REML gives as varaince estimate
simulate_IBD_blocks = 1;
plot_flag = 0;
ttt = cputime;
constant_flag = 0;
figs_dir = '../../common_disease_model/figs/IBD_sharing';
data_dir = '../../common_disease_model/data/IBD_sharing';

PLINK_file = fullfile(data_dir, 'PLINK_data');

f_vec = mat2vec(repmat([0.15 0.31 0.21 0.21 0.5 0.17 0.2 0.13 0.31 0.41 0.2 0.42 0.2 ...
    0.17 0.15 0.15 0.4 0.5 0.2 0.2 0.22 0.2 0.17 0.2], 2, 1));

%f_vec = [0.1 0.5]';
%f_vec = repmat(0.2, 10, 1); % [0.1 0.1 0.2 0.2 0.3 0.3 0.4 0.4 0.5 0.5]';

f_vec(:)=0.5; % set everything at 0.5
f_vec = repmat(f_vec, 50, 1); % take few snps to make everything fast and simple
n = length(f_vec);


if(do_one_family) % just plot one family to see that IBD-sharing estimation works
    [family_tree family_genotype_vec] = simulate_family_genotypes(1, ... % first assume diploid snps
        f_vec, iters, 'sampling', 'binary');
    family_genotype_vec = reshape(family_genotype_vec, length(f_vec), length(family_tree));
    %family_genotype_vec = reshape(family_genotype_vec, length(family_tree), length(f_vec))';
    [true_IBD_mat true_allele_sharing_mat] = kinship_coefficient(family_tree);
    
    switch snp_type
        case 'diploid'
        case 'binary' % simplified assumption - binary SNPs
            
            
    end
    
    
    m=length(family_tree)
    
    figure; imagesc(family_genotype_vec);
    figure; graph_draw(family_tree);
    figure; imagesc(true_IBD_mat); colorbar; title('True IBD-mat');
    
    
    IBD_mat = estimate_IBD_sharing_from_snps(family_genotype_vec, f_vec);
    figure; imagesc(IBD_mat); colorbar; title('Reconstructed IBD-mat');
    
    %family_genotype_vec = [0 1 1 0; 0 1 0 1]'; m=2; n=4
    family_genotype_vec = family_genotype_vec(1:2:end,:)+family_genotype_vec(2:2:end,:);
    genotype_similarity_mat = 1-(repmat(sum(family_genotype_vec), m, 1) + ...
        repmat(sum(family_genotype_vec), m, 1)' - ...
        2.*(family_genotype_vec' * family_genotype_vec))./n;
    figure; imagesc(genotype_similarity_mat); colorbar; title('Genotype similarity mat');
    
    figure; plot(true_IBD_mat(:), IBD_mat(:), '.'); % show reconstruction results
    xlabel('True IBD'); ylabel('Reconstructed IBD');
end % if do one family

if(do_slope_analysis) % here simulate many pairs of individuals with different IBD sharing
    iters = 10000; % number of pairs
    IBD_vec = randn(iters,1)*IBD_std+IBD_mean;
    IBD_vec = max(0, min(1, IBD_vec));
    
    h_all = 0.9; % heritability
    pair_genotype_vec = zeros(n,2);
    phenotype_vec = zeros(iters,2);
    estimated_IBD_vec = IBD_vec;
    estimated_IBD_const_vec = IBD_vec;
    estimated_IBD_eric_vec = IBD_vec;
    estimated_IBD_eric3_vec = IBD_vec;
    estimated_IBD_visscher_vec = IBD_vec;
    switch scenario
        case 1
            beta_vec = zeros(n,1); % all genotypes have the same influence
            beta_vec(1:n/2) = 1; % only first half of genotypes contribute
        case 2
            
    end
    
    
    
    
    for scenario=4:4 % 1:4 % loop on differenct scenarios
        switch scenario
            case 1
                
            case 2
                
            case 3
        end
        
        for i=1:iters
            if(mod(i,100)==0)
                run_iter = i
            end
            w = rand(n,1) < IBD_vec(i); % determine the set of SNPs shared IBD
            pair_genotype_vec(:,1) = rand(n,1)<f_vec;
            pair_genotype_vec(:,2) = w.*pair_genotype_vec(:,1)+...
                (1-w).*(rand(n,1)<f_vec); % simulate pair of correlated phenotypes
            for j=1:2
                phenotype_vec(i,j) = sum(pair_genotype_vec(:,j).*beta_vec);
            end
            switch scenario
                case 1 % assume we know the IBD precisely
                    estimated_IBD_vec(i) = IBD_vec(i);
                case 2 % use ALL genotypes to estimate IBD
                    tmp_estimated_mat_const = estimate_IBD_sharing_from_snps( ...
                        pair_genotype_vec, f_vec, 3+0*constant_flag); % solve for IBD_vec using ALL SNPs
                    tmp_estimated_mat = estimate_IBD_sharing_from_snps( ...
                        pair_genotype_vec, f_vec, 0+constant_flag); % solve for IBD_vec using ALL SNPs
                    tmp_estimated_mat_eric = estimate_IBD_sharing_from_snps( ...
                        pair_genotype_vec, f_vec, 2); % solve for IBD_vec using ALL SNPs
                    tmp_estimated_mat_eric3 = estimate_IBD_sharing_from_snps( ...
                        pair_genotype_vec, f_vec, 4); % solve for IBD_vec using ALL SNPs
                    tmp_estimated_mat_visscher = estimate_IBD_sharing_from_snps( ...
                        pair_genotype_vec, f_vec, 5); % solve for IBD_vec using ALL SNPs
                    
                    estimated_IBD_vec(i) = tmp_estimated_mat(1,2);
                    estimated_IBD_const_vec(i) = tmp_estimated_mat_const(1,2);
                    estimated_IBD_eric_vec(i) = tmp_estimated_mat_eric(1,2);
                    estimated_IBD_eric3_vec(i) = tmp_estimated_mat_eric3(1,2);
                    estimated_IBD_visscher_vec(i) = tmp_estimated_mat_visscher(1,2);
                    
                case 3 % use only causal SNPs to estimate IBD
                    tmp_estimated_mat = estimate_IBD_sharing_from_snps( ...
                        pair_genotype_vec(1:n/2,:), f_vec(1:n/2), constant_flag); % solve for IBD_vec using ALL SNPs
                    estimated_IBD_vec(i) = tmp_estimated_mat(1,2);
                case 4 % use only non-causal SNPs to estimate IBD
                    tmp_estimated_mat = estimate_IBD_sharing_from_snps( ...
                        pair_genotype_vec(n/2+1:end,:), f_vec(1:n/2), 2+0*constant_flag); % solve for IBD_vec using ALL SNPs
                    estimated_IBD_vec(i) = tmp_estimated_mat(1,2);
            end
        end
        n_plot = min(10000, iters);
        mean_Z = mean(phenotype_vec(:));
        sigma_Z = 0.5*(std(phenotype_vec(:,1))+std(phenotype_vec(:,2)));
        phenotype_vec = (phenotype_vec-mean_Z) ./ sigma_Z;
        phenotype_vec = phenotype_vec .* sqrt(h_all) + randn(iters,2)*sqrt((1-h_all));
        
        %    r_vec = (phenotype_vec(:,1)- phenotype_vec(:,2)).^2;
        r_vec = (phenotype_vec(:,1).*phenotype_vec(:,2));
        figure; hold on; % subplot(2,1,1); hold on;
        plot(estimated_IBD_vec(1:n_plot), r_vec(1:n_plot) , '.');
        [h_slope h_slope_stats] = regress(r_vec, [estimated_IBD_vec ones(iters,1)]);
        x_vec = linspace(min(estimated_IBD_vec), max(IBD_vec), 100);
        plot(x_vec, x_vec.*h_slope(1)+h_slope(2), 'r');
        title(['IBD sharing vs. phenotype similarity, h_{slope}^2=' num2str(h_slope(1),3) ...
            ' C.I. ( ' num2str(h_slope_stats(1,1),3) ', ' num2str(h_slope_stats(1,2), 3) ...
            '), true h_{all}^2=' num2str(h_all) ', scenario ' num2str(scenario)]);
        xlabel('IBD'); ylabel('Z_1 \times Z_2'); % \Delta(Z)');
        my_saveas(gcf, fullfile(figs_dir, ['IBD_vs_phenotype' num2str(scenario)]), ...
            {'epsc', 'fig', 'jpg', 'png'});
        
        
        figure; hold on; % subplot(2,1,2); hold on;
        plot(IBD_vec(1:n_plot), estimated_IBD_vec(1:n_plot) , '.');
        xlabel('IBD (true)'); ylabel('IBD (estiamted)');
        title(['Estiamted vs. true IBD Scenario ' num2str(scenario)]);
        if(scenario == 2)
            plot(IBD_vec(1:n_plot), estimated_IBD_const_vec(1:n_plot) , 'g.');
            plot(IBD_vec(1:n_plot), estimated_IBD_eric_vec(1:n_plot) , 'r.');
            plot(IBD_vec(1:n_plot), estimated_IBD_eric3_vec(1:n_plot) , 'm.');
            plot(IBD_vec(1:n_plot), estimated_IBD_visscher_vec(1:n_plot) , 'c.');
            legend({'MLE', 'linear (weighted)', ...
                'linear (simple average)', 'linear3', 'visscher'},4);
        end
        plot(IBD_vec(1:n_plot), IBD_vec(1:n_plot) , 'k', 'linewidth', 3);
        my_saveas(gcf, fullfile(figs_dir, ['IBD_estimators' num2str(scenario)]), ...
            {'epsc', 'fig', 'jpg', 'png'});
        
        mean((IBD_vec-estimated_IBD_vec).^2)
        mean((IBD_vec-estimated_IBD_const_vec).^2)
        mean((IBD_vec-estimated_IBD_eric_vec).^2)
        mean((IBD_vec-estimated_IBD_eric3_vec).^2)
    end % loop on scenarios
end % do slope analysis


if(do_all_pairs) % simulate a set of individuals with valid IBD matrix
    for scenario = 1:2 % 2 % try both scenarios
        if(~exist('num_individuals', 'var') || isempty(num_individuals))
            num_individuals = 10;
        end
        if(~exist('h_all', 'var') || isempty(h_all))
            h_all = 0.5; % True heritability
        end
        diag_inds = sub2ind([num_individuals num_individuals], 1:num_individuals, 1:num_individuals);
        num_iters = 2000; % repeat drawing new individuals.
        %    IBD_mat = randn(num_individuals).*IBD_std + IBD_mean;
        all_IBD_vec = []; all_phenotype_corr_vec = []; all_phenotype_vec = [];
        for i=1:num_iters % run several iterations
            run_iter = i
            old_IBD_mat=0;
            if(old_IBD_mat)
                IBD_mat = rand(num_individuals).*IBD_mean; % uniform distribution
                IBD_mat = triu(IBD_mat)+triu(IBD_mat)';
                IBD_mat = IBD_mat - diag(diag(IBD_mat)) + eye(num_individuals);
                
                % New: Create two blocks
                IBD_mat = eye(num_individuals);
                IBD_mat(1:num_individuals/2, 1:num_individuals/2) = 1-0.999925;
                IBD_mat(num_individuals/2+1:end, num_individuals/2+1:end) = 1-0.99999975;
            else % compute IBD_mat by different levels
                % New: create many levels of IBD sharing
                num_levels = 10; block_size = num_individuals/num_levels;
                IBD_mat = eye(num_individuals);
                for j=1:num_levels
                    IBD_mat((j-1)*block_size+1:end,(j-1)*block_size+1:end) = ...
                        (j-0.5)/num_levels;
                end
                IBD_mat = IBD_mat - diag(diag(IBD_mat)) + eye(num_individuals);
                IBD_u = unique(IBD_mat(:));
            end
            if(plot_flag)
                figure; imagesc(IBD_mat); colorbar;
            end
            
            %%%%%%%%%%%%%%%%  Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 z = simulate_correlated_gaussians(1000,num_individuals, IBD_mat);
            %                 z_corr_mat = z'*z;
            %                 figure; imagesc(z_corr_mat); colorbar;
            %                 if(length(IBD_u) < 500)
            %                     for i=1:length(IBD_u)
            %                         mean_z_corr(i) = mean(z_corr_mat(IBD_mat(:) == IBD_u(i)));
            %                     end
            %                 end
            %                 figure; hold on;
            %                 plot(IBD_u, mean_z_corr, 'oc', 'markersize', 15);
            %                 plot(IBD_u, mean_z_corr, '*c', 'markersize', 15);
            %%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            SNP_mat = simulate_IBD_genotypes(IBD_mat, f_vec, 1);
            snp_ids = []; people_ids = [];
            num_snps = size(SNP_mat, 1) % f_vec)
            switch scenario
                case 1
                    beta_vec = ones(num_snps,1); % all genotypes contribute equally
                case 2
                    beta_vec = ones(num_snps,1); % all genotypes have the same influence
                    beta_vec(1:num_snps/2) = 0; % only first half of genotypes contribute
            end
            num_snps = size(SNP_mat, 1) % f_vec)
            switch scenario % take only first half
                case 2
                    num_snps = num_snps/2;
            end
            %     phenotype_vec = phenotype_vec - mean(phenotype_vec);
            %     phenotype_vec = sqrt(h_all) * phenotype_vec ./ std(phenotype_vec);
            
            phenotype_vec = (SNP_mat-0.5)' * beta_vec; %  zeros(num_individuals,1);
            phenotype_vec = phenotype_vec .* 2 ./ sqrt(num_snps);
            phenotype_vec = sqrt(h_all) * phenotype_vec + randn(num_individuals, 1).*sqrt(1-h_all);
            
            all_IBD_vec = [all_IBD_vec vec2row(mat2vec(IBD_mat))];
            all_phenotype_corr_vec = ...
                [all_phenotype_corr_vec vec2row(mat2vec(phenotype_vec*phenotype_vec'))];
            all_phenotype_vec = [all_phenotype_vec vec2row(phenotype_vec)];
            
            ttt_simulate = cputime - ttt
            save_genotypes_in_PLINK_format(SNP_mat(1:num_snps,:), snp_ids, ...
                people_ids, f_vec(1:num_snps), phenotype_vec, ...
                [PLINK_file '_h_' num2str(h_all) '_beta_' num2str(scenario) ...
                '_n_' num2str(num_individuals) '_iter_' num2str(i)]); % save as plink
        end % loop on iters
        
        %        return;
        
        estimate_in_matlab = 1;
        if(estimate_in_matlab)
            ttt_estimate = cputime;
            estimated_IBD_mat = estimate_IBD_sharing_from_snps(SNP_mat, f_vec, constant_flag)
            ttt_estimate = cputime - ttt_estimate
            IBD_vec = mat2vec(IBD_mat); IBD_vec = IBD_vec(setdiff(1:num_individuals.^2, diag_inds));
            estimated_IBD_vec = mat2vec(estimated_IBD_mat); estimated_IBD_vec = estimated_IBD_vec(setdiff(1:num_individuals.^2, diag_inds));
            figure; hold on; plot(IBD_vec(:), estimated_IBD_vec(:), '.');
            
            r = corr(IBD_vec(:), estimated_IBD_vec(:));
            title(['IBD estimation (corr=' num2str(r) ')']); xlabel('true IBD'); ylabel('Estimated IBD');
            beta = polyfit(IBD_mat(:), estimated_IBD_mat(:), 1);
            x_vec = linspace(0, max(max(IBD_vec(:)), max(estimated_IBD_vec(:))), 100);
            plot(x_vec, beta(1).*x_vec + beta(2), 'r');
            plot(x_vec, x_vec, 'k');
            
            % % %             phenotype_corr_mat = phenotype_vec * phenotype_vec';
            % % %             figure; plot(IBD_mat(:), phenotype_corr_mat(:), '.'); hold on;
            % % %             [h_slope h_slope_stats] = regress(phenotype_corr_mat(:), [IBD_mat(:) ones(num_individuals^2,1)]);
            % % %             x_vec = linspace(min(IBD_mat(:)), max(IBD_mat(:)), 100);
            figure; plot(all_IBD_vec(:), all_phenotype_corr_vec(:), '.'); hold on;
            [h_slope h_slope_stats] = regress(vec2column(all_phenotype_corr_vec), ...
                [vec2column(all_IBD_vec) ones(length(all_IBD_vec),1)]);
            x_vec = linspace(min(all_IBD_vec), max(all_IBD_vec), 100);
            
            plot(x_vec, x_vec.*h_slope(1)+h_slope(2), 'r');
            title(['IBD sharing vs. phenotype similarity, h_{slope}^2=' num2str(h_slope(1),3) ...
                ' C.I. ( ' num2str(h_slope_stats(1,1),3) ', ' num2str(h_slope_stats(1,2), 3) ...
                '), true h_{all}^2=' num2str(h_all)]);
            xlabel('IBD'); ylabel('Z_1 \times Z_2'); % \Delta(Z)');
            %             mean(mat2vec(phenotype_corr_mat(1:num_individuals/2,1:num_individuals/2)))
            %             mean(mat2vec(phenotype_corr_mat(num_individuals/2+1:end,num_individuals/2+1:end)))
            %             mean(mat2vec(phenotype_corr_mat(num_individuals/2+1:end,1:num_individuals)))
            %             mean(mat2vec(phenotype_corr_mat(1:num_individuals,num_individuals/2+1:end)))
            if(length(IBD_u) < 500)
                for i=1:length(IBD_u)
                    mean_phen_corr(i) = mean(all_phenotype_corr_vec(all_IBD_vec(:) == IBD_u(i)));
                    %                    mean_phen_corr(i) = mean(phenotype_corr_mat(IBD_mat(:) == IBD_u(i)));
                end
            end
            plot(IBD_u, mean_phen_corr, '*k', 'markersize', 10);
            plot(IBD_u, mean_phen_corr, 'ok', 'markersize', 10);
            ylim([min(mean_phen_corr)-0.1, max(mean_phen_corr)+0.1]);
            
            % % %             % Hold off estimated IBD sharing
            % % %             figure; plot(estimated_IBD_mat(:), phenotype_corr_mat(:), '.'); hold on;
            % % %             [h_slope h_slope_stats] = regress(phenotype_corr_mat(:), [estimated_IBD_mat(:) ones(num_individuals^2,1)]);
            % % %             x_vec = linspace(min(IBD_mat(:)), max(IBD_mat(:)), 100);
            % % %             plot(x_vec, x_vec.*h_slope(1)+h_slope(2), 'r');
            % % %             title(['IBD sharing (estimated) vs. phenotype similarity, h_{slope}^2=' num2str(h_slope(1),3) ...
            % % %                 ' C.I. ( ' num2str(h_slope_stats(1,1),3) ', ' num2str(h_slope_stats(1,2), 3) ...
            % % %                 '), true h_{all}^2=' num2str(h_all)]);
            % % %             xlabel('IBD'); ylabel('Z_1 \times Z_2'); % \Delta(Z)');
        end % estimate in matlab
    end % loop on scenarios
end % do all pairs

if(analyze_REML_results)
    F = {'all_beta_1_h_0.5.Fhsq', 'all_beta_2_h_0.5.hsq', ...
        'all_beta_1_h_0.001.hsq', 'all_beta_2_h_0.001.hsq'};
    for i=1:length(F)
        R = loadcellfile(fullfile(data_dir, F{i}));
        h_inds = strmatch('V(1)/Vp', R(:,1));
        h_estimate(i) = mean(cell2mat(R(h_inds,2)))
        Ve_inds = strmatch('V(e)', R(:,1));
        Ve_estimate(i) = mean(cell2mat(R(Ve_inds,2)))
        Vp_inds = strmatch('Vp', R(:,1));
        Vp_estimate(i) = mean(cell2mat(R(Vp_inds,2)))
    end
    
    
end % analyze REML results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(simulate_IBD_blocks) % New: simulate regression by IBD sharing
    if(~exist('trait_type', 'var') || isempty(trait_type))
        trait_type = 'quantitative'; % 'disease'
    end
    mu = 0.25; % set disease prevalence
    mu_cc = 0.5; % set disease 'prevalence' in case-control study 
    x_mu = norminv(1-mu); % set threshold for disease
    if(machine == PC) % plot only in PC 
        plot_flag = 1;
    else
        plot_flag = 0;
    end
    
    if(~exist('num_snps', 'var') || isempty(num_snps))
        num_snps = 1000;
    end
    if(~exist('num_people_vec', 'var') || isempty(num_people_vec))
        num_people_vec = 1000; % 500:500:1000; % 1000; % :100:1000;
    end
    if(~exist('num_founders', 'var') || isempty(num_founders))
        num_founders = 10; % determines average IBD sharing
    end
    if(~exist('IBD_iters', 'var') || isempty(iters))
        iters = 100; % number of times to simulated data
    else
        iters = IBD_iters;
    end
    if(~exist('num_generations', 'var') || isempty(num_generations))
        num_generations = 10;  % determine variations in IBD sharing
    end
    num_bins = 100; % number of IBD sharing bins
    f_vec = repmat(0.5, num_snps, 1);
    num_people = num_people_vec(end); %     num_people = 2000;
    if(~exist('LP_model_str', 'var') || isempty(LP_model_str))
        LP_model_str = 'P*'; % 'debug'; % 'P*'; % 'debug'; % 'P*'; % 'debug';
    end
    switch LP_model_str
        case 'linear'
            h_pathway = 0.8; num_pathways = 1; % set LP model parameters
            h_shared_env = 0;
        case 'debug'
            h_pathway = 1; num_pathways = 10; % set LP model parameters
            h_shared_env = 0;
        case 'P*' % This one appears in the paper
            h_pathway = 0.5; num_pathways = 4; % set LP model parameters from the paper's example
            h_shared_env = 0.5 * (1-h_pathway);
        case 'P*_full' % like the one appearing in the paper but with full heritability
            h_pathway = 1; num_pathways = 4; % set LP model parameters from the paper's example
            h_shared_env = 0.5 * (1-h_pathway);
    end
    beta_snps = sqrt(num_pathways / sum(f_vec .* (1-f_vec))); % set coefficients such that variance sums to one
    model_save_str = ['LP_k_' num2str(num_pathways) '_h_' num2str(h_pathway*100) ...
        '_n_people_' num2str(num_people) '_n_snps_' num2str(num_snps) ...
        '_IBD_k0_' num2str(1/num_founders,2)];
    
    %    h_all = 1; % heritability_scale_change_MLT(h_pathway, 1, num_pathways);

    [~, ~, ~, ~, ~, h_all, tmp_h_pop] = ...
        compute_k_of_N_gaussian_statistics(0, 1, h_pathway, h_shared_env, [], ...
        'MAX', [], num_pathways, 10, 'numeric', 0, {'ACE'}); % Compute h_all (narrow sense heritability)
    h_pop = tmp_h_pop{1};
    switch trait_type % compute also heritability on disease scale
        case 'disease'
            mu_l = 1 - (1-mu)^(1/num_pathways);
            [lambda_R_numeric stat_struct] = ...
                compute_k_of_N_liabilities_statistics(num_pathways, 1, mu_l, h_pathway, h_shared_env, ...
                1, {'ACE'})
        these_are_not_the_same = h_all - stat_struct.h_liab_loci        
    end
    

    h_shared_env=0; % from now on, h_shared_env is zero for the slope analysis
    
    IBD_mean = 1/num_founders;
    IBD_res = 0.02;
    IBD_range = linspace(0, 1, 1+1/IBD_res);
    IBD_range = (IBD_range - IBD_mean) ./ (1-IBD_mean);
    
    IBD_data_file = fullfile(figs_dir, ...
        [model_save_str '_num_gen_' num2str(num_generations) '_' trait_type '.mat']);
    if(exist(IBD_data_file, 'file'))
        load_IBD_flag = 0; % temp (should be one)
    else
        load_IBD_flag = 0;
    end

    if(exist(IBD_data_file, 'file'))
        load(IBD_data_file); % load anyway so that we don't need to compute analytics again
    end
    if(load_IBD_flag)  % load everything pre-computed
    else % recompute analytics and IBD matrix
        if(~exist('rho_R_analytic', 'var') || isempty(rho_R_analytic))
            [BETA_analytic, k_sharing_vec_analytic, qtl_R_analytic, lambda_R_analytic, ~,...     % First compute everything analytically
                h_loci_analytic, h_pop_analytic] = ...
                simulate_visscher_analysis_LP(10, num_pathways, h_pathway, [], ...
                [], [], 'full-range', 'numeric', trait_type, mu, IBD_range); % make the MZ correlation the same
            switch trait_type
                case 'disease'
                    rho_R_analytic = vec2row(mu.*(vec2row(lambda_R_analytic)-1) ./ (1-mu));
                case 'quantitative'
                    rho_R_analytic = vec2row(qtl_R_analytic);                    
            end
            IBD_sharing_vec_analytic = k_sharing_vec_analytic .* (1-IBD_mean) + IBD_mean; % transfer from k-IBD to blocks-IBD
            save(IBD_data_file, 'BETA_analytic', 'k_sharing_vec_analytic', 'rho_R_analytic', ...     % First compute everything analytically
                'h_loci_analytic', 'h_pop_analytic', 'IBD_sharing_vec_analytic');
        end
        IBD_sharing_vec_analytic = k_sharing_vec_analytic .* (1-IBD_mean) + IBD_mean; % transfer from k-IBD to blocks-IBD
        
        [input_IBD_mat estimated_IBD_mat ...  % Simulate one IBD mat and many SNP mats
            SNP_mat founder_SNP_mat input_SNP_founder_identity_mat ...
            input_founder_identities input_poiss_points] = ...
            simulate_IBD_blocks_sharing(num_founders, num_generations, num_people, f_vec);
        
        save(IBD_data_file, 'BETA_analytic', 'k_sharing_vec_analytic', 'rho_R_analytic', ... %  'qtl_R_analytic', ...     % First compute everything analytically
            'h_loci_analytic', 'h_pop_analytic', 'input_IBD_mat', 'estimated_IBD_mat', ...  % Simulate one IBD mat and many SNP mats
            'SNP_mat', 'founder_SNP_mat', 'input_SNP_founder_identity_mat', ...
            'input_founder_identities', 'input_poiss_points', 'IBD_sharing_vec_analytic');
    end % if load IBD
    
    all_LP_phen_vec_I = []; all_LP_phen_vec_J = []; 
    all_LP_phen_vec_I_normalized = []; all_LP_phen_vec_J_normalized = []; 
    all_SNP_mat_I = []; all_SNP_mat_J = [];
    
    [snp_start_inds snp_end_inds] = ...
        divide_region_to_blocks(1, num_snps, ceil(num_snps/num_pathways));
    
    beta_phenotype_IBD_vec = [];
    beta_phenotype_IBD_range_vec_std_theoretical = [];
    beta_phenotype_IBD_vec_local = [];
    beta_phenotype_IBD_vec_quad = [];
    IBD_all_bins_vec = 0:0.01:1; % set bins
    all_iters_smoothed_phenotype_corr_vec{1} = zeros(num_bins+1,1);
    all_iters_smoothed_phenotype_corr_vec_counts{1} = zeros(num_bins+1,1);
    all_iters_smoothed_phenotype_corr_vec_std{1} = zeros(num_bins+1,1); % compute std in each bin
    
    num_reg = IBD_mean * num_bins; % different slopes in different intervals
    bin_slope = zeros(num_reg,iters);
    IBD_a = 0.49; IBD_b = 0.51;
    [I_IBD_bin J_IBD_bin] = find((input_IBD_mat < IBD_b) & (input_IBD_mat > IBD_a));
    mean_prevalence=0; % record mean prevalence 
    for i_iter=1:iters % loop on iters
        sprintf('Run iter=%ld out of %ld', i_iter, iters)
        time_simulate_snps=cputime;
        [IBD_mat estimated_IBD_mat ...  % Simulate one IBD mat and many SNP mats
            SNP_mat founder_SNP_mat SNP_founder_identity_mat] = ...
            simulate_IBD_blocks_sharing(num_founders, num_generations, num_people, f_vec, ...
            input_IBD_mat, input_founder_identities, input_poiss_points, input_SNP_founder_identity_mat);
        SNP_mat = SNP_mat - repmat(f_vec, 1, num_people); % normalize genotypes to have mean zero
        IBD_mat = IBD_mat - diag(diag(IBD_mat)); % eye(num_people);
        estimated_IBD_mat = estimated_IBD_mat - diag(diag(estimated_IBD_mat)); eye(num_people);
        time_simulate_snps = cputime - time_simulate_snps
        
        time_simulate_phenotype = cputime;
        P = randperm(num_snps); phen_SNP_mat = SNP_mat(P,:); % permute not to have adjacent SNPs in the same pathway !!! 
        LP_phenotype_vec = zeros(num_people,1)-99999999999; % start with very negative
        for i=1:num_pathways
            temp_pathay_vec = ...
                sum(beta_snps .* phen_SNP_mat(snp_start_inds(i):snp_end_inds(i),:))' .* sqrt(h_pathway) + ...
                randn(num_people,1) .* sqrt(1-h_pathway);% Generate phenotype
            LP_phenotype_vec = max(LP_phenotype_vec, temp_pathay_vec);
        end
        mu_vec(i_iter) = mean(LP_phenotype_vec);
        std_vec(i_iter) = std(LP_phenotype_vec);
        all_LP_phen_vec_I = [all_LP_phen_vec_I vec2row(LP_phenotype_vec(I_IBD_bin))];
        all_LP_phen_vec_J = [all_LP_phen_vec_J vec2row(LP_phenotype_vec(J_IBD_bin))];
        LP_phenotype_vec = LP_phenotype_vec - mean(LP_phenotype_vec); % normalize
        LP_phenotype_vec = LP_phenotype_vec ./ std(LP_phenotype_vec);

        switch trait_type % convert to disease scale (AFTER normalization!)
            case 'disease'
                LP_phenotype_vec = double(LP_phenotype_vec > x_mu); % simulate disease according to liability threshold model
        end
        mean_prevalence = mean_prevalence + mean(LP_phenotype_vec);
        
        
        all_LP_phen_vec_I_normalized = [all_LP_phen_vec_I_normalized vec2row(LP_phenotype_vec(I_IBD_bin))];
        all_LP_phen_vec_J_normalized = [all_LP_phen_vec_J_normalized vec2row(LP_phenotype_vec(J_IBD_bin))];
% % %         if(length(all_SNP_mat_I) < 2000) % avoid memory problems
% % %             all_SNP_mat_I = [all_SNP_mat_I phen_SNP_mat(:, I_IBD_bin)];
% % %             all_SNP_mat_J = [all_SNP_mat_J phen_SNP_mat(:, J_IBD_bin)];
% % %         end
        
        
        LP_phenotype_corr_mat = LP_phenotype_vec * LP_phenotype_vec';
        LP_phenotype_corr_mat = LP_phenotype_corr_mat - diag(diag(LP_phenotype_corr_mat)); % remove diagonal

% % %         LP_phenotype_disease_corr_mat = LP_phenotype_disease_vec * LP_phenotype_disease_vec'; % compute the same for disease 
% % %         LP_phenotype_disease_corr_mat = LP_phenotype_disease_corr_mat - diag(diag(LP_phenotype_disease_corr_mat)); % remove diagonal
        
        time_simulate_phenotype = cputime - time_simulate_phenotype
        
        smooth_span = 1000; res = smooth_span/2;
        IBD_sharing_bins = 100;
        
        %    IBD_std = 1/num_generations; % That's not quite true. Compute empirically !
        IBD_std = std(squareform(IBD_mat));
        beta_phenotype_IBD_vec_std_theoretical = 1 ./ (sqrt(num_people_vec.^2./2)*IBD_std);
        
        
        time_loop_on_people = cputime;
        for i_num_people = 1:length(num_people_vec)
            cur_num_people = num_people_vec(i_num_people);
            res = max(1, round(cur_num_people^2/10000));
            model_str = ['LP(k=' num2str(num_pathways), ', h_{pathway}^2=' num2str(100*h_pathway) '%).' ...
                ' n_{people}=' num2str(cur_num_people) ' n_{causal-snps}=' num2str(num_snps) ...
                ', IBD: k_0=' num2str(1/num_founders) ', \sigma(k)=' num2str(IBD_std,2)  '. h_{all}^2=' num2str(h_all,3)];
            
            %    IBD_vec = IBD_vec + randn(1,length(IBD_vec)).*0.0000001; % make them distinct
            for j=1:1 % loop on regressing on true vs. estimated IBD vec
                if(j == 1) % take true IBD mat
                    IBD_vec{j} = squareform(IBD_mat(1:cur_num_people,1:cur_num_people)); % must be positive !!!
                    IBD_str = 'True IBD';
                else % take estimated IBD mat (could be negative!!!)
                    IBD_vec{j} = squareform(estimated_IBD_mat(1:cur_num_people,1:cur_num_people));
                    IBD_str = 'Estimated IBD';
                end % if j==1
                [IBD_vec{j} IBD_sort_perm] = sort(IBD_vec{j}); % sort for smoothing data
                LP_phenotype_corr_vec{j} = ... % copy all correlation
                    squareform(LP_phenotype_corr_mat(1:cur_num_people,1:cur_num_people));

                switch trait_type % normalize correlation. But we need to do this for each run!!!!
                    case 'disease' % correct: rho_R = (lambda_R - mu^2)/ (mu*(1-mu))
                        cur_mu = mean(LP_phenotype_vec);
                        LP_phenotype_corr_vec{j} = ... % correct correlation
                            (LP_phenotype_corr_vec{j}-cur_mu^2) ./ (cur_mu*(1-cur_mu)); % This just gives the correlation coefficient
                end                
                LP_phenotype_corr_vec{j} = LP_phenotype_corr_vec{j}(IBD_sort_perm); % sort by IBD sharing 
                
                [beta_phenotype_IBD{j} beta_phenotype_IBD_interval{j}] = ...
                    regress(LP_phenotype_corr_vec{j}', ...
                    [ones(length(IBD_vec{j}), 1) IBD_vec{j}']); % fit a linear regression
                [beta_linear{j}] = polyfit(IBD_vec{j}', LP_phenotype_corr_vec{j}', 1); % linear regression again
                IBD_std = std(squareform(IBD_mat));
                %                IBD_range = [max(0, IBD_mean - 0.5*IBD_std), min(1, IBD_mean + 0.5*IBD_std)];

                
                
                for i_range = 1:3 % try different ranges
                    IBD_estimation_range(i_range,:) = [IBD_mean-IBD_mean.*i_range/3, IBD_mean+IBD_mean.*i_range/3];
                    f_range = find((IBD_vec{j} < IBD_estimation_range(i_range,2)) & ...
                        (IBD_vec{j} >= IBD_estimation_range(i_range,1))); % exclude zero! % include zero !
                    [beta_linear_near_IBD{j}(i_range,:)] = polyfit(IBD_vec{j}(f_range)', ...
                        LP_phenotype_corr_vec{j}(f_range)', 1); % local linear regression
                    %                    [beta_quad{j}(i_range,:)] = polyfit(IBD_vec{j}(f_range)', ...
                    %                        LP_phenotype_corr_vec{j}(f_range)', 2); % local quadratic regression
                    num_points_in_range{j}(i_range,i_iter) = length(f_range); % see how many points did we loose
                    IBD_std_in_range{j}(i_range)  = std(IBD_vec{j}(f_range));
                end
                
                x_vec = min(IBD_vec{j}):0.01:1;
                
                if(j == 1) % change this. We run only j==1 
                    smoothed_phenotype_corr_vec{j} = zeros(num_bins+1,1);
                    smoothed_phenotype_corr_vec_std{j} = zeros(num_bins+1,1);
                    smoothed_phenotype_corr_vec_counts{j} = zeros(num_bins+1,1);
                    cur_bins_vec = 1:max(floor(IBD_vec{j}'*num_bins)+1);
                    smoothed_phenotype_corr_vec{j}(cur_bins_vec) = ...
                        accumarray(floor(IBD_vec{j}'*num_bins)+1, LP_phenotype_corr_vec{j}', [], @mean);
                    smoothed_phenotype_corr_vec_std{j}(cur_bins_vec) = ...
                        accumarray(floor(IBD_vec{j}'*num_bins)+1, LP_phenotype_corr_vec{j}', [], @std);
                    smoothed_phenotype_corr_vec_counts{j}(cur_bins_vec) = ...
                        accumarray(floor(IBD_vec{j}'*num_bins)+1, LP_phenotype_corr_vec{j}', [], @length);
                    IBD_bins_vec = [0:num_bins] ./ num_bins;
                    
                    
                    for i_reg = 1:num_reg % regress on close lines
                        [~, bin_middle] = min(abs(IBD_bins_vec - IBD_mean));
                        bin_inds = [bin_middle-i_reg:bin_middle+i_reg];
                        bin_beta = regress( ...
                            smoothed_phenotype_corr_vec{1}(bin_inds), ...
                            [IBD_bins_vec(bin_inds)' ones(length(bin_inds), 1)]);
                        bin_slope(i_reg,i_iter) = bin_beta(1);
                    end
                    
                    
                    
                    
                    %
                    %                     [zeros(1, num_bins) ...
                    %                         linspace(min(IBD_vec{j}), max(IBD_vec{j}), ...
                    %                         length(smoothed_phenotype_corr_vec{j})-num_bins)];
                else % use estimated IBD mat - can be negative
                    smoothed_phenotype_corr_vec{j} = ...
                        accumarray(floor(IBD_vec{j}'*num_bins)+num_bins, ...
                        LP_phenotype_corr_vec{j}', [], @mean);
                    smoothed_phenotype_corr_vec_std{j} = ...
                        accumarray(floor(IBD_vec{j}'*num_bins)+num_bins, ...
                        LP_phenotype_corr_vec{j}', [], @std);
                    IBD_bins_vec = [zeros(1, num_bins) ...
                        linspace(min(IBD_vec{j}), max(IBD_vec{j}), ...
                        length(smoothed_phenotype_corr_vec{j})-num_bins)];
                end
                smoothed_phenotype_corr_vec2{j} = smooth(LP_phenotype_corr_vec{j}, smooth_span);  % Smooth data to see a trend
                
                all_iters_smoothed_phenotype_corr_vec{j} = ...
                    all_iters_smoothed_phenotype_corr_vec{j} + ...
                    smoothed_phenotype_corr_vec{j} .* smoothed_phenotype_corr_vec_counts{j}; % Compute cumulative total
                all_iters_smoothed_phenotype_corr_vec_std{j} = ...
                    all_iters_smoothed_phenotype_corr_vec_std{j} + ...
                    smoothed_phenotype_corr_vec_std{j}.^2 .* smoothed_phenotype_corr_vec_counts{j} + ...
                    smoothed_phenotype_corr_vec{j}.^2 .* smoothed_phenotype_corr_vec_counts{j}; % compute variance
                
                all_iters_smoothed_phenotype_corr_vec_counts{j} = ...
                    all_iters_smoothed_phenotype_corr_vec_counts{j} + ...
                    smoothed_phenotype_corr_vec_counts{j}; % Compute cumulative counts
                
                if((i_iter <= 2) && (i_num_people == length(num_people_vec)))% plot estimation for one iteration
                    %                    figure; plot(squareform(IBD_mat), squareform(estimated_IBD_mat), '.')
                    %                    [I J] = find(IBD_mat > 0.99)
                    % % %                     figure; hold on; % plot(estimated_IBD_mat(:), LP_phenotype_corr_mat(:), '.');
                    % % %                     xlabel('Estimated IBD mat'); ylabel('Phenotypic correlation mat');
                    % % %                     plot(x_vec, x_vec .* beta_phenotype_IBD{j}(2) + beta_phenotype_IBD{j}(1), 'r');
                    % % %                     plot(IBD_vec{j}, smoothed_phenotype_corr_vec2{j}, '.g');
                    if((i_iter==1) && (i_num_people==length(num_people_vec))) % if num people
                        if(plot_flag)
                            figure; hold on;
                            hist(squareform(IBD_mat), 100);
                            xlabel('IBD'); ylabel('Freq.');
                            title(['IBD sharing distribution. \mu_{IBD}=' ...
                                num2str(mean(squareform(IBD_mat)),2) ', \sigma_{IBD}=' ...
                                num2str(std(squareform(IBD_mat)),2)]);
                        end % if plot
                    end
                    
                    if(plot_flag)
                        figure; hold on;
                        xlabel([IBD_str ' mat']); ylabel('Phenotypic correlation mat');
                        plot(IBD_vec{j}(1:res:end), smoothed_phenotype_corr_vec2{j}(1:res:end), '.m');
                        plot(IBD_sharing_vec_analytic, rho_R_analytic, 'k', 'linewidth', 3);
                        
                        plot(IBD_bins_vec, smoothed_phenotype_corr_vec{j}, 'b', 'linewidth', 2);
                        %                     errorbar(IBD_bins_vec, ...
                        %                         smoothed_phenotype_corr_vec{j}, smoothed_phenotype_corr_vec_std{j}, ...
                        %                         'linewidth', 2); % , '.g');
                        plot(x_vec, x_vec .* beta_phenotype_IBD{j}(2) + beta_phenotype_IBD{j}(1), ...
                            'r', 'linewidth', 2);
                        %                     for i_range = 1:3
                        %                         plot(x_vec, x_vec.^2 .* beta_quad{j}(i_range, 1) + x_vec .* beta_quad{j}(i_range, 2) + ...
                        %                             beta_quad{j}(i_range,3), 'g', ...
                        %                             'linewidth', 3, 'linestyle', symbol_vec{i_range});
                        %                     end
                        for i_range = 1:3
                            plot(x_vec, x_vec .* beta_linear_near_IBD{j}(i_range,1) + ...
                                beta_linear_near_IBD{j}(i_range,2), 'g', ... % green clearer
                                'linewidth', 3, 'linestyle', symbol_vec{i_range});
                        end
                        % '], \beta (quad.)=' num2str(beta_quad{j}(2),3) ...
                        title(['IBD vs. \rho(Z_1,Z_2). ' model_str ...
                            ', \sigma_{\beta} (theor.)=' ...
                            num2str(beta_phenotype_IBD_vec_std_theoretical(end),3) ...
                            ', \beta = ' num2str(beta_phenotype_IBD{j}(2),3) ' [' ...
                            num2str(beta_phenotype_IBD_interval{j}(2,1),3) ', ' ...
                            num2str(beta_phenotype_IBD_interval{j}(2,2),3) ...
                            '], \beta (range)=' num2str(beta_linear_near_IBD{j}(1),3)]);
                        xlim([0,1]); ylim([-0.1,1.1]);
                        % %                         ['local-quad-fit (\beta=' num2str(beta_quad{j}(1,2),3) ')'], ...
                        % %                         ['local-quad-fit (\beta=' num2str(beta_quad{j}(2,2),3) ')'], ...
                        % %                         ['local-quad-fit (\beta=' num2str(beta_quad{j}(3,2),3) ')'], ...
                        legend({'smoothed-data', 'theoretical-\rho',  'bin-average-error', ...
                            ['linear-fit (\beta=' num2str(beta_phenotype_IBD{j}(2),3) ')'], ...
                            ['local-linear-fit [' num2str(IBD_estimation_range(1,1),2) ', ' num2str(IBD_estimation_range(1,2),2) ...
                            '] (\beta=' num2str(beta_linear_near_IBD{j}(1,1),3) ')'], ...
                            ['local-linear-fit [' num2str(IBD_estimation_range(2,1),2) ', ' num2str(IBD_estimation_range(2,2),2) ...
                            '] (\beta=' num2str(beta_linear_near_IBD{j}(2,1),3) ')'], ...
                            ['local-linear-fit [' num2str(IBD_estimation_range(3,1),2) ', ' num2str(IBD_estimation_range(3,2),2) ...
                            '] (\beta=' num2str(beta_linear_near_IBD{j}(3,1),3) ')']}, 2);
                        if(i_iter == 1) % save only one example
                            my_saveas(gcf, fullfile(figs_dir, ...
                                ['IBD_sharing_vs_phenotypic_correlation_model_' model_save_str]), ...
                                {'epsc', 'fig', 'jpg', 'png'});
                        end
                    end % if plot
                end % if iter=1
                beta_phenotype_IBD_vec{j}(i_iter, i_num_people) = beta_phenotype_IBD{j}(2); % collect all betas
                for i_range=1:3
                    %                    beta_phenotype_IBD_vec_quad{j,i_range}(i_iter, i_num_people) = ...
                    %                        2.*beta_quad{j}(i_range,1) .* IBD_mean + beta_quad{j}(i_range,2); % compute slope at k0
                    beta_phenotype_IBD_vec_local{j,i_range}(i_iter, i_num_people) = ...
                        beta_linear_near_IBD{j}(i_range,1);
                    beta_phenotype_IBD_range_vec_std_theoretical{j}(i_range,:) = ...
                        1 ./ (sqrt(num_people_vec.^2./2)*IBD_std_in_range{j}(i_range));  %need to take into account # points available
                end
            end % loop on j (true vs. estimated IBD vec)
        end % loop on number of people
        time_loop_on_people = cputime - time_loop_on_people
    end % loop on iterations
    mean_prevalence = mean_prevalence/iters
    
    % Plot average over all iterations
    all_iters_smoothed_phenotype_corr_vec{1} = ...
        all_iters_smoothed_phenotype_corr_vec{1} ./ ...
        all_iters_smoothed_phenotype_corr_vec_counts{1};
    all_iters_smoothed_phenotype_corr_vec_std{1} = ...
        sqrt( (all_iters_smoothed_phenotype_corr_vec_std{1} - ...
        all_iters_smoothed_phenotype_corr_vec{1}.^2) ./ ...
        all_iters_smoothed_phenotype_corr_vec_counts{1} );
    all_iters_smoothed_phenotype_corr_vec_std{1} = ...
        all_iters_smoothed_phenotype_corr_vec_std{1} ./ ...
        sqrt(all_iters_smoothed_phenotype_corr_vec_counts{1});
    
    
% % % %     switch trait_type % normalize correlation. But we need to do this for each run!!!! 
% % % %         case 'disease' % correct: rho_R = (lambda_R - mu^2)/ (mu*(1-mu))
% % % %             all_iters_smoothed_phenotype_corr_vec{1} = ... % correct correlation 
% % % %                 (all_iters_smoothed_phenotype_corr_vec{1}-mu^2) ./ (mu*(1-mu)); % This just gives the correlation coefficient           
% % % % %            all_iters_smoothed_phenotype_corr_vec_std{1} = ...
% % % % %                (all_iters_smoothed_phenotype_corr_vec_std{1}-mu^2) ./ (mu*(1-mu));  % st.d. is scaled the same way?? why? 
% % % %     end
    
    average_file_name = fullfile(figs_dir, ...
        ['IBD_sharing_vs_phenotypic_correlation_model_' model_save_str '_iters_' num2str(iters) '_average']);
    if(plot_flag) % plot results of averaging
        figure; hold on; 
        %title(['IBD vs. \rho(Z_1,Z_2). ' model_str ' iters=' num2str(iters) ' average']);
        plot(IBD_sharing_vec_analytic, rho_R_analytic, 'k', 'linewidth', 3);
        plot_beta_ind = min(5,size(bin_slope,1)); % which interval to take 

        % 
        plot_beta =  mean(bin_slope(plot_beta_ind,:));
        
        rho_correction_factor = 1; % mu*(1-mu) / normpdf(x_mu)^2; % assume population sampling

        switch trait_type % normalize slope and multiply by constant to get an estimator for heritability
            case 'quantitative'
                plot_mean_h_slope = plot_beta*(1-IBD_mean);
                plot_beta_std = std(bin_slope(plot_beta_ind,:))*(1-IBD_mean);
            case 'disease'
                plot_mean_h_slope = plot_beta*(1-IBD_mean) * ...
                    rho_correction_factor; % assume population sampling
                plot_beta_std = std(bin_slope(plot_beta_ind,:))*(1-IBD_mean) * ...
                    rho_correction_factor; % assume population sampling. Std. is scaled the same way! 
        end
        
        plot_y0 = -plot_beta .* IBD_mean; % make sure regression line passes at (k0, 0)
        good_bins = find(all_iters_smoothed_phenotype_corr_vec_counts{1} >= 1000)
        errorbar(IBD_bins_vec(good_bins), all_iters_smoothed_phenotype_corr_vec{1}(good_bins), ...
            all_iters_smoothed_phenotype_corr_vec_std{1}(good_bins), ...
            'b', 'linewidth', 2);
        plot(x_vec, plot_beta .* x_vec + plot_y0, 'r', 'linewidth', 2); % plot regression line
        
        %            1 ./ sqrt(all_iters_smoothed_phenotype_corr_vec_counts{1}), ...
        %            'b', 'linewidth', 2);
        xlabel('IBD sharing'); ylabel('\rho(IBD)');
        legend('\rho (Gauss. approx.)', '\rho empirical', 'slope at k_0', 2);           legend boxoff;
        y_lim = get(gca, 'ylim')
        line([IBD_mean IBD_mean], [y_lim(1)  0], ...
            'color', 'k', 'linestyle', '--');
        line([0 IBD_mean], [0 0], ...
            'color', 'k', 'linestyle', '--');
        text(IBD_mean, y_lim(1)*1.1, 'k_0', 'fontsize', 14); 
        xlim([0,1]);
        switch trait_type
            case 'quantitative'
                %                 slope_str = ['$\beta(1-k_0) =' num2str(plot_mean_h_slope,3) ...
                %                     ' \pm ' num2str(plot_beta_std,3) '$'];
                slope_str = ['$h_{slope(mean-IBD)}^2 =' num2str(plot_mean_h_slope,3) ...
                    ' \pm ' num2str(plot_beta_std,3) '$'];
                stat_struct.h01_loci2 = -1;
            case 'disease'
                %                 slope_str = ['$\frac{\beta(1-k_0) \mu^2(1-\mu)^2}' ...
                %                     '{\varphi(x_{\mu})^2 \mu_{cc}(1-\mu_{xx})} =' num2str(plot_mean_h_slope,3) ...
                %                     ' \pm ' num2str(plot_beta_std,3) '$'];
                slope_str = ['$h_{slope(mean-IBD)}^2 =' num2str(plot_mean_h_slope,3) ...
                    ' \pm ' num2str(plot_beta_std,3) '$'];
        end
        text(IBD_mean*2, IBD_mean*0, slope_str, 'fontweight', 'bold', 'Interpreter','latex');
        
    title([model_str '. iters=' num2str(iters) ' average, h_{all}^2=' ...
        num2str(h_all,3) ', h_{01}^2=' num2str(stat_struct.h01_loci2,3)]);
%    title([repmat(' ', 1, 80) 'supp-fig6'], 'fontsize', 16, 'fontweight', 'bold');  
        my_saveas(gcf, average_file_name, {'epsc', 'fig', 'jpg', 'png'});        
    
        full_figure; hold on; % New figure; plot mean and st.d. as a function of interval size
        beta_bin_mean = mean(bin_slope,2)*(1-IBD_mean);
        beta_bin_std = std(bin_slope,[], 2)*(1-IBD_mean);
        interval_size_vec = 2*diff(IBD_bins_vec(1:2)) .* (1:length(beta_bin_mean));
        errorbar(interval_size_vec, beta_bin_mean, beta_bin_std, 'linewidth', 2);
        plot(interval_size_vec, repmat(h_all, length(interval_size_vec), 1), 'k--', 'linewidth', 2);
        legend('estimator', 'true h_{all}^2');         legend boxoff;
        xlabel('interval size'); ylabel('(1-k_0) \beta estimator'); 
        title([model_str '. iters=' num2str(iters) ' average']);
        xlim([0 max(interval_size_vec)+0.01]);
        my_saveas(gcf, [average_file_name '_different_intervals'], {'epsc', 'fig', 'jpg', 'png'});
        
    end % if plot
    
    num_reg = IBD_mean * num_bins; % different slopes in different intervals
    R_slope = cell(num_reg,7); bin_beta = zeros(num_reg,2);
    for i_reg = 1:num_reg % regress on close lines
        [~, bin_middle] = min(abs(IBD_bins_vec - IBD_mean))
        bin_inds = [bin_middle-i_reg:bin_middle+i_reg];
        bin_beta(i_reg,:) = regress(all_iters_smoothed_phenotype_corr_vec{1}(bin_inds), ...
            [IBD_bins_vec(bin_inds)' ones(length(bin_inds), 1)]);
        R_slope{i_reg,1} = ['[' num2str(IBD_bins_vec(bin_inds(1)),2) ...
            '-' num2str(IBD_bins_vec(bin_inds(end)),2) ']'];
        R_slope{i_reg,2} = num2str( ...
        100*sum(all_iters_smoothed_phenotype_corr_vec_counts{1}(bin_inds)) / ...
            (iters*num_people*(num_people-1)/2), 3);        
        R_slope{i_reg,3} = num2str(bin_beta(i_reg,1),4);
        R_slope{i_reg,4} = num2str(bin_beta(i_reg,1)*(1-IBD_mean),4);
        
        R_slope{i_reg,5} = num2str(mean(bin_slope(i_reg,:)),4); % mean of many regression
        R_slope{i_reg,6} = num2str(mean(bin_slope(i_reg,:))*(1-IBD_mean),4); % std of many regression
        R_slope{i_reg,8} = num2str(std(bin_slope(i_reg,:))*(1-IBD_mean),4); % std of many regression
    end
    R_slope = [{'interval', '% points', 'combined \beta (slope)', 'combined  \beta*(1-k0) (h_all estimator)', ...
        '', '\beta mean', '\beta*(1-k0) mean', '\beta*(1-k0) st.d.'}' R_slope']';
    R_slope = [['Regression:', repmat({''}, 1, 7)]' R_slope']';
    
    
    
    R_ave = [IBD_bins_vec' all_iters_smoothed_phenotype_corr_vec_counts{1} ...
        all_iters_smoothed_phenotype_corr_vec{1} ...
        all_iters_smoothed_phenotype_corr_vec_std{1}];
    R_ave = [{'bin', 'num-points', 'mean', 'std'}' num2str_cell(num2cell(R_ave), 7)']';
    
    R_slope = [R_slope' cell(size(R_ave,1) - size(R_slope,1), size(R_slope,2))']';
    R_ave = [R_ave R_slope];
    savecellfile(R_ave, [average_file_name '.txt']);
    
    beta_phenotype_IBD_vec_mean = cell(1,1);
    beta_phenotype_IBD_vec_std = cell(1,1);
    beta_phenotype_IBD_vec_local_mean = cell(1,1);
    beta_phenotype_IBD_vec_local_std = cell(1,1);
    for j=1:1 % Plot the average beta with it's variation
        beta_phenotype_IBD_vec_mean{j} = mean(beta_phenotype_IBD_vec{j},1) .* (1-IBD_mean);
        beta_phenotype_IBD_vec_std{j} = std(beta_phenotype_IBD_vec{j},[],1) .* (1-IBD_mean);
        for i_range=1:3
            %             beta_phenotype_IBD_vec_quad_mean{j}(i_range,:) = mean(beta_phenotype_IBD_vec_quad{j,i_range},1) .* (1-IBD_mean);
            %             beta_phenotype_IBD_vec_quad_std{j}(i_range,:) = std(beta_phenotype_IBD_vec_quad{j,i_range},1) .* (1-IBD_mean);
            beta_phenotype_IBD_vec_local_mean{j}(i_range,:) = mean(beta_phenotype_IBD_vec_local{j,i_range},1) .* (1-IBD_mean);
            beta_phenotype_IBD_vec_local_std{j}(i_range,:) = std(beta_phenotype_IBD_vec_local{j,i_range},1) .* (1-IBD_mean);
        end
        frac_points_in_range_vec{j} = mean(num_points_in_range{j},2) ./ (num_people*(num_people-1)/2);
    end
    
    
    x_slope_vec = IBD_sharing_vec_analytic(1:end-1) + 0.5*diff(IBD_sharing_vec_analytic);
    y_slope_vec = (1-IBD_mean) .* diff(vec2row(rho_R_analytic))./diff(IBD_sharing_vec_analytic);
    y_high_res_slope_vec = interp(y_slope_vec, 100);
    x_high_res_slope_vec = interp(x_slope_vec, 100);
    
    if(plot_flag)
        figure; hold on; % plot slope as function of IBD sharing
        plot(x_slope_vec, y_slope_vec, 'linewidth', 2);
        plot(IBD_sharing_vec_analytic, repmat(h_all, length(IBD_sharing_vec_analytic), 1), ...
            'k', 'linewidth', 2);
        line([IBD_mean IBD_mean],  [0 h_all], 'linestyle', '--', 'color', 'k', 'linewidth', 2);
        title('Slope as function of IBD');
        xlabel('IBD'); ylabel('(1-k_0) d \rho(IBD) / d IBD');
        legend('slope', 'h_{all}^2');
        my_saveas(gcf, fullfile(figs_dir, ...
            ['IBD_sharing_exact_slope_model_' model_save_str]), ...
            {'epsc', 'fig', 'jpg', 'png'});
        
        
        figure; hold on; % Plot beta estimation
        plot(num_people_vec, repmat(h_all, length(num_people_vec), 1), 'k-', 'linewidth', 5);
        errorbar(num_people_vec, ...
            beta_phenotype_IBD_vec_mean{1}, beta_phenotype_IBD_vec_std{1}, 'linewidth', 2);
        
        %     for i_range=1:3 % plot quad
        %         errorbar(num_people_vec+i_range*10, ...
        %             beta_phenotype_IBD_vec_quad_mean{1}(i_range,:), beta_phenotype_IBD_vec_quad_std{1}(i_range,:), ...
        %             'r', 'linewidth', 2, 'linestyle', symbol_vec{i_range});
        %     end
        for i_range=1:3 % plot empirical estimators
            errorbar(num_people_vec+i_range*20, ...
                beta_phenotype_IBD_vec_local_mean{1}(i_range,:), beta_phenotype_IBD_vec_local_std{1}(i_range,:), ...
                'c', 'linewidth', 2, 'linestyle', symbol_vec{i_range});
        end
        %     errorbar(num_people_vec, ...
        %         repmat(h_all, length(num_people_vec), 1), beta_phenotype_IBD_vec_std_theoretical, ...
        %         'g', 'linewidth', 2);
        for i_range=1:3 % plot theoretical estimators
            errorbar(num_people_vec+i_range*20-10, ...
                repmat(h_all, length(num_people_vec), 1), beta_phenotype_IBD_range_vec_std_theoretical{j}(i_range,:), ...
                'g', 'linewidth', 2, 'linestyle', symbol_vec{i_range});
        end
        
        
        %    errorbar(num_people_vec, ... % for now drop estimated IBD (assume we know it precisely)
        %        beta_phenotype_IBD_vec_mean{2}, beta_phenotype_IBD_vec_std{2}, 'r');
        
        title(['h_{slope}^2. ' model_str '. iters=' num2str(iters)]);
        xlabel('num people'); ylabel('\beta estimator');
        % %         ['\beta-local-quad1 [' num2str(IBD_mean*(1-1/3),2) '..' num2str(IBD_mean*(1+1/3),2) '], #=' ...
        % %         num2str(frac_points_in_range_vec{1}(1)*100,2) '%'], ...
        % %         ['\beta-local-quad2 [' num2str(IBD_mean*(1-2/3),2) '..' num2str(IBD_mean*(1+2/3),2) '], #=' ...
        % %         num2str(frac_points_in_range_vec{1}(2)*100,2) '%'], ...
        % %         ['\beta-local-quad3 [' num2str(IBD_mean*(1-3/3),2) '..' num2str(IBD_mean*(1+3/3),2) '], #=' ...
        % %         num2str(frac_points_in_range_vec{1}(3)*100,2) '%'], ...
        legend_vec = {'h_{all}^2', '\beta (IBD) (\pm st.d.)'};
        for i_range=1:3
            legend_vec = [legend_vec ...
                ['\beta-local [' num2str(IBD_mean*(1-i_range/3),2) '..' ...
                num2str(IBD_mean*(1+i_range/3),2) '], #=' ...
                num2str(frac_points_in_range_vec{1}(i_range)*100,2) '%']]; %, ...
        end
        for i_range=1:3
            legend_vec = [legend_vec ...
                ['\beta theoretical (\pm st.d.) [' num2str(IBD_mean*(1-i_range/3),2) '..' ...
                num2str(IBD_mean*(1+i_range/3),2) '], #=' ...
                num2str(frac_points_in_range_vec{1}(i_range)*100,2) '%']];
        end
        %     ['\beta-local2 [' num2str(IBD_mean*(1-2/3),2) '..' num2str(IBD_mean*(1+2/3),2) '], #=' ...
        %         num2str(frac_points_in_range_vec{1}(2)*100,2) '%'], ...
        %         ['\beta-local3 [' num2str(IBD_mean*(1-3/3),2) '..' num2str(IBD_mean*(1+3/3),2) '], #=' ...
        %         num2str(frac_points_in_range_vec{1}(3)*100,2) '%'], ...
        %         '\beta theoretical (\pm st.d.)'}); % '\beta (estimated-IBD) (\pm st.d.)',
        legend(legend_vec);
        ylim([-0.1,2]);
        my_saveas(gcf, fullfile(figs_dir, ...
            ['IBD_sharing_h_slope_estimator_model_' model_save_str]), ...
            {'epsc', 'fig', 'jpg', 'png'});
        
        figure; hold on; plot(x_slope_vec, y_slope_vec, 'linewidth', 3);
        plot(x_high_res_slope_vec, y_high_res_slope_vec, 'r');
    end % if plot
    
    
    [~, I_sibs] = min(abs(x_high_res_slope_vec - 0.5));
    h_sib = y_high_res_slope_vec(I_sibs) .* (1-IBD_mean);
    
    IBD_sharing_vec_analytic_high_res = interp(IBD_sharing_vec_analytic, 100);
    rho_R_analytic_high_res = interp(vec2row(rho_R_analytic), 100);
    [~, I_left] = min(abs(IBD_sharing_vec_analytic_high_res - IBD_mean));
    
    beta_phenotype_IBD_vec_mean_expected = ...
        (1-IBD_mean) * (rho_R_analytic_high_res(end) - rho_R_analytic_high_res(I_left)) / ...
        (IBD_sharing_vec_analytic_high_res(end) - IBD_sharing_vec_analytic_high_res(I_left));
    for i_range=1:3
        [~, I_left] = min(abs(x_high_res_slope_vec - IBD_estimation_range(i_range,1)));
        [~, I_right] = min(abs(x_high_res_slope_vec - IBD_estimation_range(i_range,2)));
        beta_phenotype_IBD_vec_local_mean_expected(i_range) = ...
            mean(y_high_res_slope_vec([I_left I_right]));
    end
    precision =3;
    R = strsplit(' ', strdiff(strdiff(model_str, '{'), '}'));
    R = [[R{1,1} R{1,2}] R(1,3:end) ...
        ['h_pop^2=' num2str(h_pop,precision)] ['h_sib-slope^2=' num2str(h_sib,precision)] ...
        ['iters=' num2str(iters)]]; % Take model parameters
    R{2,1} = 'estimator\#individuals';
    R{3,2} = 'Range'; R{3,3} = '% points'; R{3,4} = '\sigma(IBD)';
    for j=1:length(num_people_vec)
        R{2,j*5+1} = num2str(num_people_vec(j));
        R{3,j*5+1} = 'Exp.'; R{3,j*5+2} = 'Obs.';
        R{3,j*5+3} = '(+/-) std. (exp.)'; R{3,j*5+4} = '(+/-) std. (obs.)';
    end
    
    for j=1:length(num_people_vec) % Copy the global estimator
        R{4,1} = 'linear'; % estimator name
        R{4,2} = '[0,1]';
        R{4,3} = '100%';
        R{4,4} = num2str(IBD_std,2);
        R{4,j*5+1} = num2str(beta_phenotype_IBD_vec_mean_expected,precision); % expected
        R{4,j*5+2} = num2str(beta_phenotype_IBD_vec_mean{1}(j),precision); % observed
        R{4,j*5+3} = num2str(beta_phenotype_IBD_vec_std_theoretical(j),precision); % +/- st.d. (expected)
        R{4,j*5+4} = num2str(beta_phenotype_IBD_vec_std{1}(j),precision); % +/- st.d. (observed)
    end
    for i_range=1:3 % copy the different local estimators
        R{i_range+4,1} = 'local-linear'; % estimator name
        R{i_range+4,2} = ['[' num2str(IBD_estimation_range(i_range,1),precision) ', ' ...
            num2str(IBD_estimation_range(i_range,2),precision) ']'];
        R{i_range+4,3} = [num2str(100*frac_points_in_range_vec{1}(i_range),precision) '%'];
        R{i_range+4,4} = num2str(IBD_std_in_range{1}(i_range),precision);
        for j=1:length(num_people_vec)
            R{i_range+4,j*5+1} = num2str(beta_phenotype_IBD_vec_local_mean_expected(i_range),precision); % expected slope. What is it?
            R{i_range+4,j*5+2} = num2str(beta_phenotype_IBD_vec_local_mean{1}(i_range,j),precision); % observed
            R{i_range+4,j*5+3} = num2str(beta_phenotype_IBD_vec_local_std{1}(i_range,j),precision); % +/- st.d. (expected)
            R{i_range+4,j*5+4} = num2str(beta_phenotype_IBD_range_vec_std_theoretical{1}(i_range,j),precision); % +/- st.d. (observed)
        end
    end
    R{end+1,1} = ''; R{end+1,1} = ''; % add spaces
    
    savecellfile(R, fullfile(figs_dir, ...     % Create also a table and save
        ['IBD_sharing_h_slope_estimator_model_' model_save_str '.txt']));
    
    
    % %    figure; imagesc(IBD_mat-diag(diag(IBD_mat))); colorbar; figure; isposdef(IBD_mat)
    % %     figure; hold on; plot(IBD_mat(:), estimated_IBD_mat(:), '.');
    % %     [C C_pval] = corr(IBD_mat(:), estimated_IBD_mat(:));
    % %     beta = regress(estimated_IBD_mat(:), [ones(length(IBD_mat(:)), 1) IBD_mat(:)]);
    % %     x_vec = 0:0.01:1;
    % %     plot(x_vec, x_vec .* beta(2) + beta(1), 'r');
    % %     title(['IBD estimation corr=' num2str(C,3) ', n_{SNPs}=' num2str(num_snps)]);
    % %     xlabel('True IBD'); ylabel('Estimated IBD');
    % %     figure; imagesc(founder_SNP_mat); colorbar;  title('Founder SNPs');
    % %     figure; subplot(2,1,1); imagesc(SNP_founder_identity_mat'); colorbar; title('Founder Inheritance Regions');
    % %     subplot(2,1,2);  imagesc(SNP_mat'); colorbar; title('Population SNPs');
    
end % simulate IBD blocks

ttt = cputime - ttt

ttt2 = ttt;
create_ppt=0;
if(create_ppt)% Create a powerpoint presentation
    fig_files = GetFileNames(fullfile(figs_dir, '*.jpg'), 1)
    mxdom2ppt(fig_files,'ppt')
    publish(fig_files{1}, 'ppt')
end



num_points = 20000; K=10; h_pathway=4/9; c = zeros(K,1);
x = randn(K, num_points) .*sqrt(h_pathway);
x_env1 = randn(K, num_points).*sqrt(1-h_pathway);
x_env2 = randn(K, num_points).*sqrt(1-h_pathway);
y1 = max(x+x_env1); y2 = max(x+x_env2);
for i=1:K
    c(i) = corr(x(i,:)'+x_env1(i,:)', y1');
end
h_is = sum(c.^2)  .* h_pathway
rho_mz_twins = corr(y1', y2')






% Simulate data to see the correlations of LP phenotypes 
h_pathway = 0.5; % 0.99999;
num_snps_per_pathway = 100; % Now simulate everything discrete
x_d = binornd(round(num_snps_per_pathway*h_pathway), 0.5, K, num_points); % common genetic part
x_env1_d = binornd(round(num_snps_per_pathway*(1-h_pathway)), 0.5, K, num_points);
x_env2_d = binornd(round(num_snps_per_pathway*(1-h_pathway)), 0.5, K, num_points); % unique environmental part
x_d = x_d - 0.5*round(num_snps_per_pathway*h_pathway);
x_env1_d = x_env1_d - 0.5*round(num_snps_per_pathway*(1-h_pathway));
x_env2_d = x_env2_d - 0.5*round(num_snps_per_pathway*(1-h_pathway));


y1_g = sum(x_d+x_env1_d); y2_g = sum(x_d+x_env2_d);
y1_d = max(x_d+x_env1_d); y2_d = max(x_d+x_env2_d); 
c_d = zeros(K,1);
for i=1:K
    c_d(i) = corr(x_d(i,:)'+x_env1_d(i,:)', y1_d');
end
h_d_is = sum(c_d.^2)  .* h_pathway
rho_dz_twins_d = corr(y1_d', y2_d')

corr_in_snps = corr(y1_g', y2_g')



% % % % Debug LP(10). These should be 0.26!!! Not 0.4 !!! 
% % % LP_corrs_normalized = corr(all_LP_phen_vec_I_normalized', all_LP_phen_vec_J_normalized')
% % % LP_corrs = corr(all_LP_phen_vec_I', all_LP_phen_vec_J')
% % % 
% % % num_pairs = size(all_SNP_mat_I,2)
% % % 
% % % for i=1:num_pairs
% % %     if(mod(i,100)==0)
% % %         run_corr_i = i
% % %     end
% % %     all_LP_SNP_corr(i) = corr(all_SNP_mat_I(:,i), all_SNP_mat_J(:,i)); 
% % % end
% % % 
% % % 
% % % all_LP_phen_vec_I_again = zeros(1,num_pairs)-999999999999;
% % % all_LP_phen_vec_J_again = all_LP_phen_vec_I_again;
% % % for i=1:10
% % %     x_dd_I(i,:) = sum(all_SNP_mat_I((i-1)*100+1:i*100,:));
% % %     x_dd_J(i,:) = sum(all_SNP_mat_J((i-1)*100+1:i*100,:));
% % %     
% % %     all_LP_phen_vec_I_again = max(all_LP_phen_vec_I_again, ...
% % %         sum(all_SNP_mat_I((i-1)*100+1:i*100,:)));
% % %     all_LP_phen_vec_J_again = max(all_LP_phen_vec_J_again, ...
% % %         sum(all_SNP_mat_J((i-1)*100+1:i*100,:)));
% % % end
% % % y_dd_I = sum(x_dd_I); y_dd_J = sum(x_dd_J);
% % % y_dd_phen_I = max(x_dd_I); y_dd_phen_J = max(x_dd_J);
% % % corr(all_LP_phen_vec_I_again', all_LP_phen_vec_J_again')
% % % corr(y_dd_phen_I', y_dd_phen_J')
% % % figure; hist(all_LP_SNP_corr, 100);
% % % mean(all_LP_SNP_corr)
% % % 
% % % 
% % % figure; hist(y1_d, min(y1_d):max(y1_d))
% % % figure; hist(y2_d, min(y2_d):max(y2_d))
% % % 
% % % figure; hist(y_dd_phen_I, min(y_dd_phen_I):max(y_dd_phen_I))
% % % figure; hist(y_dd_phen_J, min(y_dd_phen_J):max(y_dd_phen_J))
% % % 
% % % figure; hist2d_draw([y1_d' y2_d'], -5:25, -5:25, [], [], 'y1_d sim.');
% % % figure; hist2d_draw([y_dd_phen_I' y_dd_phen_J'], -5:25, -5:25, [], [], 'LP sim.');
