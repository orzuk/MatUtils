% Compute power parameters that should appear in figure and in supp. info.
% of the paper.
% There are numerous tests we run. We focus on marginal, pairwise and one
% pathway test
%
% Input:
% Why are there so many different h values? 
% h_x - total heritability explained in each pathway 
% h_x_one_locus_vec - heritability of one SNP on the disease scale ??? (assuming LT)
%
function  compute_heritability_power_parameters(isoheritability_flag, ...
    power_pathway_alpha_vec, ...
    power_iters, n_samples_vec, ...
    h_x, h_x_one_locus_vec,  ...
    freq, K, N, mu, mu_ind, mu_vec, num_loci, ...
    test_name_vec, test_type_vec, test_stat_vec, sampling_type_vec, try_test_inds)


%[p_x_y p_x_z p_x_x_z disease_grr_vec] = ...
%    heritability_to_p_z_x_MLT(cur_h_x_one_locus_LT, freq, K, N, mu);
% power_marginal(:,mu_ind)  = ... % Compute power
%     compute_association_power(mat2vec(p_x_z') , ...
%     n_samples_vec, [], power_alpha, ... % power_alpha, ... % Should be power_alpha !!
%     [], ... % we can afford more iters in marginal (faster)
%     'armitage', 'chi-square-analytic', 'case-control');
%         [power_pairwise_epistasis p_vals_epistasis stats_epistasis] = ... % compute epistasis power: takes a long time ...
%             compute_association_power(p_x_x_z, ...
%             n_samples_vec, [], power_pairwise_alpha, ...
%             power_iters, ... % use less iterations (this is slower)
%             'epistasis', 'probit', 'case-control'); % run PAIRWISE epistasis (should be low power)


% Compute pathway-based power. This is more complicated. First, for
% each n_sample compute how many true loci were detected. Then, add
% their heritability and determine their grr and interaction. Then,
% compute power.
% % % % % % % %         plot_discovered_loci_to_epistasis_iterative=0;
% % % % % % % %         if(plot_discovered_loci_to_epistasis_iterative)
% % % % % % % %             for i=length(n_samples_vec):-1:1 % go backwards so heavy memory is run first
% % % % % % % %                 run_pathway_power = i
% % % % % % % %                 [p_x_y p_x_z p_x_x_z pathway_grr_vec(i)] = heritability_to_p_z_x_MLT( ...
% % % % % % % %                     h_x_one_locus*num_loci*power_marginal(i), ...
% % % % % % % %                     freq, K, N, mu);
% % % % % % % %
% % % % % % % %                 % compute 2x2 tables for MLT model
% % % % % % % %                 [cur_power_pathway_epistasis p_vals_pathway_epistasis{i} stats_pathway_epistasis{i}] = ... % compute epistasis power: takes a long time ...
% % % % % % % %                     compute_association_power(p_x_x_z, ...
% % % % % % % %                     [min(n_samples_vec(i)-10,10000):100000:n_samples_vec(i) n_samples_vec(i)], ...
% % % % % % % %                     [], power_pathway_alpha, power_iters, ... % use less iterations (this is slower)
% % % % % % % %                     'epistasis', 'logistic', 'case-control'); % simulate blocks to avoid memory problems
% % % % % % % %                 power_pathway_epistasis(i) = cur_power_pathway_epistasis(end);
% % % % % % % %             end % loop on num sample
% % % % % % % %         end % plot epistasis of discovered loci

%    load(['power_MLT_data_iso_' num2str(isoheritability_flag)]);
%    n_samples_vec = 2*round(logspace(2,4,10)/2); %   % 2*round(logspace(2,6,50)/2); % try short n samples vec (run faster)

h_x_MLT_one_locus_vec = zeros(mu_ind, length(h_x_one_locus_vec)); 
non_centrality_parameter = cell(max(try_test_inds), mu_ind); 
for i=1:length(h_x_one_locus_vec) % loop on heritability of a single locus/pathway 
    compute_power_for_h = i
    sprintf('try_power_computation,h_x=%2.1f%% ,mu=%2.2f%%', ...
        100*h_x_one_locus_vec(i), mu*100)
    h_x_MLT_one_locus_vec(mu_ind,i) = ... % Transfer from hertiability explained to variance explained !!! 
        heritability_scale_change_MLT(N*h_x_one_locus_vec(i), 1, N, mu, 'MLT'); % .*h_x
    
    h_x_MLT_one_locus_vec(mu_ind,i)=h_x_MLT_one_locus_vec(mu_ind,i)/2;
    [p_x_y p_x_z p_x_x_z] = heritability_to_p_z_x_MLT( ...
        h_x_MLT_one_locus_vec(mu_ind,i)/1, freq, K, N, mu);
    [p_x_y_marginal p_x_z_marginal p_x_x_z_marginal] = heritability_to_p_z_x_MLT( ...
        h_x_MLT_one_locus_vec(mu_ind,i)/1, freq, K, N, mu); % use half effect size for Armitage test

    h_x_one_locus_again = p_z_x_marginal_to_heritability(mat2vec(p_x_z')') % Test that variance didn't change
    
    % Alternative: 
%    grr = heritability_to_genetic_relative_risk(h_x_one_locus_vec(i), 'liability', freq, mu);
%     p_x_z_alternative = vec2mat(genetic_relative_risk_to_p_z_x_marginal(freq, grr, mu), 2)';
    p_full_sample_mat = [0.5 sqrt(h_x_MLT_one_locus_vec(mu_ind,i)) N mu]; % the MLT model
    p_full_sample_mat_LT = [0.5 sqrt(h_x_one_locus_vec(i)) 1 mu]; % compare to LT model
    model_params = [h_x_one_locus_vec(i) h_x_MLT_one_locus_vec(mu_ind,i) mu N freq];
    p_pathway_vec = {p_full_sample_mat, p_full_sample_mat_LT, ...
        p_full_sample_mat, p_full_sample_mat_LT, ...
        p_full_sample_mat, p_full_sample_mat_LT, p_full_sample_mat, ...
        p_full_sample_mat, p_full_sample_mat_LT, p_full_sample_mat, ...
        mat2vec(p_x_z_marginal')', mat2vec(p_x_z_marginal')', mat2vec(p_x_z_marginal')', ...
        p_x_x_z, p_x_x_z, p_x_x_z}; % need also to insert pairwise epistasis here !
    
    for j=try_test_inds % 1:num_test_stats % run all tests. Currently heaviest is chi-square in bins
        %                run_test_stat = j
        sprintf(['try_pathway_epistasis,h_x=%2.1f%% ,mu=%2.2f%% test_stat=%ld ' test_name_vec{j}], ...
            100*h_x_one_locus_vec(i), mu*100, j)
        [power_pathway_epistasis_cell{j,mu_ind}(i,:), ~, ~, tmp_ncp]= ... % compute epistasis power: takes a long time ...
            compute_association_power(p_pathway_vec{j}, ...
            n_samples_vec, [], power_pathway_alpha_vec(j), ...
            power_iters, ... % use less iterations (this is slower)
            test_type_vec{j}, test_stat_vec{j}, sampling_type_vec{j}, [], model_params);    % simulate in blocks to avoid meory problems
        if(~isempty(tmp_ncp))
            non_centrality_parameter{j,mu_ind}(i,:) = tmp_ncp;
        end
    end
    
end % loop on different valeus of h_x

% % figure; plot(n_samples_vec, non_centrality_parameter{j,mu_ind});
% % title('NCP as function of sample size');
% % xlabel('n samples'); ylabel('NCP'); 
% % legend( num2str_cell(num2cell(100*h_x_one_locus_vec)),4);
% % 
% % figure; hold on; 
% % title('NCP as function of effect size (var. explained)');
% % plot(100*h_x_one_locus_vec, non_centrality_parameter{j,mu_ind}(:,1)) % only (:,1) means ... 
% % plot(100*h_x_one_locus_vec, non_centrality_parameter{j,mu_ind}(:,1), '*')
% % xlabel('h_x^2'); ylabel('NCP'); 



% % % % % %     AssignGeneralConstants;
% % % % % %     mu=0.01;   n_samples_vec = 2*round(logspace(2,7,50)/2); % 2*round(logspace(2,6,50)/2);
% % % % % %     N=2; h_x_one_locus_vec = logspace(-3,-1,7); % possible heritabilities on X axis
% % % % % %     pow_mat = zeros(length(h_x_one_locus_vec), length(n_samples_vec));
% % % % % %     for i=1:length(h_x_one_locus_vec)% Here we compute the power analytically (a temp hack)
% % % % % %         [gof_mu gof_sigma pow_mat(i,:)] = ...
% % % % % %             MLT_goodness_of_fit_stat_moments(h_x_one_locus_vec(i), mu, 2, 1, ...
% % % % % %         n_samples_vec, 0.01, 'MLT', 100^2);
% % % % % %         close all;
% % % % % %     end
% % % % % %     figure; semilogx(n_samples_vec, pow_mat, 'linewidth', 2);
% % % % % %     xlabel('Sample size'); ylabel('detection power (%)');
% % % % % %     legend_vec = num2str_cell(num2cell(h_x_one_locus_vec*100),2,[],1);
% % % % % %     legend(legend_vec, 2);
% % % % % %     title(['power for \mu=' num2str(mu) ' for different h_x^2. Significance \alpha = 0.01']);
% % % % % %     my_saveas(gcf, ['../../common_disease_model/figs/epistasis_power/pathway_epistasis_gof_mu_' ...
% % % % % %         num2str(mu)], format_fig_vec);

save(['power_LP_disease_' num2str(N) '_1_data_iso_' num2str(isoheritability_flag) ], ...
    'n_samples_vec', ... 'power_marginal', ...% 'power_pairwise_epistasis', ... % 'power_pathway_epistasis', ...
    'N', 'K', 'mu_vec', 'freq', 'num_loci', ... % 'pathway_grr_vec', ...
    'power_pathway_epistasis_cell', 'try_test_inds', ...
    'test_name_vec', 'test_type_vec', 'test_stat_vec', 'sampling_type_vec', ...
    'power_pathway_alpha_vec',  'non_centrality_parameter', ...
    'h_x', 'h_x_one_locus_vec', 'h_x_MLT_one_locus_vec', 'power_iters');  % save data (heavy run ..)
%     'power_alpha', 'power_pairwise_alpha', 'power_pathway_alpha', ...

