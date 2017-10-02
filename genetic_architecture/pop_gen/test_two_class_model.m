% Test different computations for the two-class model: power-calculations, likelihood
AssignGeneralConstants;
AssignRVASConstants;
total_test_time = cputime;
figs_dir = '../../common_disease_model/figs/EyreWalker';
new_figs_dir = fullfile(figs_dir, 'new_eric');
tau_vec=1; % [0.1 0.5 1 2 10];
sigma_epsilon_vec=[0 1]; % 0.0000000000000000000001; % noise in coupling between fitness and trait
delta=1; % multiplicative coefficient relating fitness and trait
theta_vec= [3000 300 30 3]; % Gamma scale parameter
k=0.2; % Gamma shape parameter. Value given by Eyre-Walker
num_points=5000000; % number of points to simulate
N = 10000; % effective population size
in_matlab_flag = 1;


global cumsum_log_vec;

test_power = 1; run_compute_power=0; % here test detection power.
test_ratio_likelihood = 0; % test rare alleles using maximum likelihood for enrichment data
test_likelihood = 0; poisson_flag = 0; % test rare alleles using maximum likelihood for full data (genotype+phenotype)
compare_power_to_var_explained = 0; % just plot power vs. variance explained 


% demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/sim_results_set2/mat';
%demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/files_var_eric/mat'; % new file !!
%demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/EuropeFixed/files_europe_fixed/files_var_eric/mat'; % new file (corrected Europe) !!
demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/EuropeFixed/correct_fine_s_high_res/mat/all'; % new file (corrected Europe) !!


% % % demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/EuropeFixed/fine_s_var/out_var/mat'; % new file (corrected Europe fine-s values) !!

demographic_models_struct = [];
%demographic_models_struct.model_str = {'european', 'varselection1', 'varselection2', 'twophase'}; % different population models
demographic_models_struct.file_names = GetFileNames(fullfile(demographic_sims_dir, 'cumul*.mat'), 1);  % NEW! add before/after bottle-neck




for i=1:length(demographic_models_struct.file_names)
    demographic_models_struct.model_str{i} = ...
        strdiff(strdiff(strdiff(remove_suffix_from_file_name(remove_dir_from_file_name(demographic_models_struct.file_names{i})), ...
        'cumul_'), '_all'), '_summary');
end

% % % {fullfile(demographic_sims_dir, 'cumul_europ.mat'), ...
% % %     fullfile(demographic_sims_dir, 'cumul_varsel1.mat'), ...
% % % fullfile(demographic_sims_dir, 'cumul_varsel1.mat'), ...
% % %     fullfile(demographic_sims_dir, 'cumul_varsel2.mat'), ...
% % %     fullfile(demographic_sims_dir, 'cumul_2phase.mat')};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(test_power) % test detection power in rare-variant GWAS
    p_val_cutoff_vec = [0.05 / 20000]; % we don't need nominal significance for now % 0.05 
    p_val_cutoff = 0.05 / 20000; % account for # of genes in the genoms
    alpha_vec = 0.3; % 1/3; % 1/2; %%% 1/3; % 0.999999999; % proportion of functional 'null' rare mutations
    beta_vec = linspace(0.001, 10.001, 201); % logspace(-6,1,100); % effect size
    beta_vec = [beta_vec -linspace(0.0001, 0.5, 201)]; % NEW !!! add negarive effects !!! 
    prevalence = 0.05; % prevalence (assume 5%)
    mu=mu_per_site; % 1.6*10^(-8); % mutation rate per nucleotide per generation
    gene_length = 625; % 1500; % typical gene length in nucleotides. Makes the gene mutation rate mu_g at 10^(-5)
    
    
    n_cases_vec = 1000:1000:100000; n_controls_vec = n_cases_vec; n_vec = n_cases_vec + n_controls_vec;
    [~, i_beta] = min(abs(beta_vec  - 3)); % 61; % choose specific beta
    [~, i_n] = min(abs(n_vec - 50000)); %  25; % choose specific n
    
    rare_cumulative_per_gene = 1; % cumulative allele frequency of rare alleles per gene
    f_rare = 0.01; % frequency below which an allele is considered 'rare'
    f_rare_vec = logspace(-6,0,1000); % different possible maximal values
    f_rare_vec(end) = 0.99999999999999;
    s_null_vec = logspace(-6,-1,500); % try different selection coefficients
    show_s_null = [0 logspace(-5, -1, 9)]; % go up to 10^-1%     [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01];
    for i=1:length(show_s_null)
        [~, show_s_null_ind(i)] = min(abs(s_null_vec - show_s_null(i)));
    end
    s_null = 0.01; % selection coefficient on null (functional) alleles - pretty strong
    
    
    full_enrichment_alpha_vec = [alpha_vec 1]; % represent how many of alleles are null
    w_x_null_mat = cell(2,1); w_x_harmless = cell(2,1); w_all = cell(2,1);
    c_cumulative = cell(2,1); frac_null_by_freq_cumulative = cell(2,1);
    equilibrium_parameters_output_file = fullfile(new_figs_dir, 'equilibrium', 'two_class_equilibrium_parameters');
    if(run_compute_power) %   > 99999)
        for i=1:2  % need to transfer this also to parallel jobs
            compute_alpha_i = i
            [two_class_stat_struct{i}, w_x_null_mat{i}, w_x_harmless{i}, w_all{i} c_cumulative{i}, frac_null_by_freq_cumulative{i}] = ...
                compute_two_class_model_parameters(s_null_vec, ...
                f_rare_vec, full_enrichment_alpha_vec(i), rare_cumulative_per_gene, N);
        end
        save(equilibrium_parameters_output_file, ...
            's_null_vec', 'f_rare', 'f_rare_vec', 'beta_vec', 'alpha_vec', 'full_enrichment_alpha_vec', ...
            'two_class_stat_struct', 'w_x_null_mat', 'w_x_harmless', 'w_all', 'c_cumulative', 'frac_null_by_freq_cumulative', ...
            'n_controls_vec', 'n_cases_vec', 'n_vec', 'rare_cumulative_per_gene', ...
            'N', 'prevalence', 'gene_length', 'mu', ...
            'p_val_cutoff', 'p_val_cutoff_vec', 'show_s_null_ind', ...
            'demographic_models_struct');
    else
        save(equilibrium_parameters_output_file, '-append', 'demographic_models_struct'); % TEMP!
    end % if run_compute_power
    load(equilibrium_parameters_output_file);
    title_str = ['Detect. power rare. \alpha=' num2str(alpha_vec*100,3) '%' ...
        ' functional, f^*=' num2str(f_rare*100,2) '%, c=' ...
        num2str(rare_cumulative_per_gene*100,3) '% cum. freq., p-val cutoff = ' num2str(p_val_cutoff, 3)];
    
    
    if(machine == PC)
        if(0 > 99999)
            save_prevalence = prevalence; prevalence = 0.01;
            plot_two_class_equilibrium_statistics(two_class_stat_struct{1}, w_x_null_mat, w_x_harmless, w_all, ... % two_class_stat_struct, w_x_null_mat, w_x_harmless, w_all,
                f_rare_vec, frac_null_by_freq_cumulative, c_cumulative, ...
                N, show_s_null, show_s_null_ind, s_null_vec, rare_cumulative_per_gene, alpha_vec,  prevalence, ...
                title_str, figs_dir, new_figs_dir, [], [], demographic_models_struct); % Plot properties of distribution - this should also be part of the figures for paper
            prevalence = save_prevalence;
        end
        in_matlab_flag = 1;
    else
        in_matlab_flag = 0;
    end
    
    %     [N_half_mat N_half_mat_only_null_enrichment N_half_mat_by_enrichment ...
    %         power_mat new_power_mat n_samples_half_power_vec f_rare_power_vec] = ...
    close all;
    if(run_compute_power)  % NEW! run each demographic model separately in order to save time !!!
        demographic_models_struct.model_str{end+1} = 'equilibrium';
        demographic_models_struct.num_models = length(demographic_models_struct.model_str);
        
        for i_d = 2:2 % 1:demographic_models_struct.num_models
            power_parameters_output_file = fullfile(new_figs_dir, demographic_models_struct.model_str{i_d}, ...
                'two_class_power_parameters.mat');
            
            %             [num_samples_struct power_struct] = ...
            %                 compute_two_class_power_parameters(equilibrium_parameters_output_file, ...
            %                 power_parameters_output_file, i_d);
            power_job_str = ['[num_samples_struct power_struct] = ' ...
                'compute_two_class_power_parameters(''' equilibrium_parameters_output_file ''', ' ...
                '''' power_parameters_output_file ''', ' num2str(i_d) ');']
            
            if(in_matlab_flag)
                eval(power_job_str);
            else % submit job
                SubmitMatlabJobToFarm(power_job_str, ['out/compute_power_' demographic_models_struct.model_str{i_d} '.out'], ...
                    'priority', [], [], 4);
            end
            
            
            % % % %             two_class_stat_struct, ...
            % % % %                 s_null_vec, f_rare, f_rare_vec, beta_vec, alpha_vec, full_enrichment_alpha_vec, ...
            % % % %                 n_controls_vec, n_cases_vec, n_vec, rare_cumulative_per_gene, ...
            % % % %                 w_x_null_mat, w_x_harmless, ...
            % % % %                 N, prevalence, gene_length, mu, ...
            % % % %                 p_val_cutoff, p_val_cutoff_vec, show_s_null_ind, c_cumulative, frac_null_by_freq_cumulative, ...
            % % % %                 power_parameters_output_file, demographic_models_struct, i_d); % compute paremters representing power
        end % loop on demographic models
    end % if run_compute_power
    
    
    
    % Need to move this into the plot function !!!!
    %    load(power_parameters_output_file); % just load instead of recomputing all power parameters, which is time-consuming
    
    new_alpha_vec = 0.01:0.005:1; % This measures how many of the rare variants are functional
    new_kept_alleles_vec = 0.01:0.01:1; % This measures how many alleles we kept in enrichment
    
    if(machine == PC)
        plot_two_class_power_statistics(equilibrium_parameters_output_file, new_figs_dir, new_figs_dir);
        % % % % % %         num_samples_struct, power_struct, ...
        % % % % % %             two_class_stat_struct, w_x_null_mat, w_x_harmless, w_all, ...
        % % % % % %             frac_null_by_freq_cumulative, ...
        % % % % % %             f_rare, f_rare_vec, beta_vec, ...
        % % % % % %             s_null_vec, show_s_null_ind,  i_beta, i_n, ...
        % % % % % %             title_str, figs_dir, new_figs_dir, ...
        % % % % % %             n_vec, new_kept_alleles_vec, alpha_vec, new_alpha_vec, ...
        % % % % % %             N, prevalence, gene_length, mu, ...
        % % % % % %             rare_cumulative_per_gene,  demographic_models_struct); % plot power curves
        %     f_rare, f_rare_vec, f_rare_power_vec, beta_vec, ...
        %             N_half_mat, N_half_mat_only_null_enrichment, ...
        %             s_null_vec, show_s_null_ind,  i_beta, i_n, ...
        %             title_str, figs_dir, new_figs_dir, n_samples_half_power_vec, power_mat, new_power_mat, N_half_mat_by_enrichment, ...
        %             n_vec, new_kept_alleles_vec, alpha_vec, new_alpha_vec, ...
        %             rare_cumulative_per_gene); % plot power curves
    end
end % if test_power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test absorption-time (two versions should agree !!!)
%absorption_time_by_selection(s, theta, N, x_min, x_max, weight_flag)
% N=10000;
% T = absorption_time_by_selection(0.0001, 1, N, 0.0001, 0.01, 'freq')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(test_ratio_likelihood)  % Here test ratio of ...???
    t1_vec = [0.25 0.5 0.95 0.05]; t2_vec = [0.75 0.5 0.95 0.05]; R_vec = [6 3]; % enrichment
    alpha_vec = 0.01:0.01:1; % fraction of deletirious mutations
    true_alpha = 0.3; true_beta = 1; % true effect size and fraction of deleterious mutations.
    n1 = 1000; n2 = 1000; % number of individuals in each tail
    f = 0.03; % allele frequency in population
    r1 = 78 ; r2 = 50; t1 = 0.25; t2 = 0.75; alpha = 0.3;
    beta_vec = NaN(length(alpha_vec), length(t1_vec)); V = beta_vec;
    s1_vec = norminv(t1_vec); s2_vec = norminv(t2_vec);
    for i=1:length(t1_vec) % compute different alpha-beta curves
        s1_vec(i) = MixtureOfGaussiansICDF(t1_vec(i), [1-true_alpha*f, true_alpha*f], [0, true_beta], [1, 1]);
        s2_vec(i) = MixtureOfGaussiansICDF(t2_vec(i), [1-true_alpha*f, true_alpha*f], [0, true_beta], [1, 1]);
        f1 = (f / t1_vec(i)) * ( true_alpha * normcdf(s1_vec(i)-true_beta) + (1-true_alpha) * normcdf(s1_vec(i)) ); % need to normalize
        f2 = (f / (1-t2_vec(i))) * ( 1 - true_alpha * normcdf(s2_vec(i)-true_beta) - (1-true_alpha) * normcdf(s2_vec(i)) ); % need to normalize
        r1 = f1*2*n1; r2 = f2*2*n2; % compute enrichment
        [~, ~, ~, alpha_min(i)] = ratio_QTL_to_var_explained(n1, n2, r1, r2, t1_vec(i), t2_vec(i), 1, f);
        good_alpha_inds = find(alpha_vec > alpha_min(i));
        
        recover_true_beta = ratio_QTL_to_var_explained(n1, n2, r1, r2, t1_vec(i), t2_vec(i), true_alpha, f);
        
        [beta_vec(good_alpha_inds,i), ~, V(good_alpha_inds,i)] = ...
            ratio_QTL_to_var_explained(n1, n2, r1, r2, t1_vec(i), t2_vec(i), alpha_vec(good_alpha_inds), f);
    end
    legend_vec = [repmat('t_1=', length(t1_vec), 1) num2str(t1_vec') repmat(', t_2=', length(t1_vec), 1)  num2str(t2_vec')];
    
    %
    figure; plot(alpha_vec, beta_vec); hold on;
    plot(true_alpha, true_beta, 'k*', 'markersize', 10);
    legend(legend_vec); xlabel('\alpha'); ylabel('\beta');
    xlim([min(alpha_min) 1]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(test_likelihood) % New: test likelihood
    L = 100; % number of loci (gene length)
    n = 2000; % number of individuals
    mu = mu_per_site; % mutation rate per-nucleotide
    trait_struct = []; 
    trait_struct.type = 'disease'; % 'disease'; % simulate either disease or quqantitative traits
    trait_struct.prevalence = 0.1; % disease prevalence for disease traits
    rare_cumulative_per_gene = 1; % cumulative allele frequency of rare alleles per gene (theta)
    target_size_by_class_vec = 10000 .* [5 5 5]; %[5 1 10]; % [synonymous,  stop, missense]
    theta = 4*N*mu; % effective mutation rate
    iters = 20; % number of simulations to perform
    full_flag = 0; % use summary statistics (faster)
    
    test_two_class_model_likelihood(N, L, n, mu, trait_struct, rare_cumulative_per_gene, target_size_by_class_vec, ...
        theta, poisson_flag, iters, full_flag, figs_dir);
    
    
end % test likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_compute_eric=0;
if(tmp_compute_eric) % Temp computation for Eric: m_s(x)
    N=10000; s_vec =logspace(-6, -2, 100); ctr=1; y_vec = []; f_vec = [1/(2*N) 0.001 0.01 0.1 1];
    for ff = f_vec
        y_vec(ctr,:) = absorption_time_by_selection(s_vec, 1, N, 0, min(ff, 0.9999999999), 'freq') * 2*N;
        ctr=ctr+1;
    end
    figure; loglog(s_vec, y_vec, 'linewidth', 2); legend(num2str(f_vec')); title('Num. alleles with freq. <= f / num. alleles at birth');
    xlabel('-s'); ylabel('m_s(f)'); my_saveas(gcf, 'proportion_of_rare_alleles_normalized_by_birth_rate', {'fig', 'pdf'});
end % tmp compute for eric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





if(compare_power_to_var_explained)
    grr_vec = 1 + logspace(-6, 3, 100);
    mu_vec = [0.1]; %  0.05 0.1]; % change prevalence
    f_vec = [0.0001 0.001  0.01 0.05 0.1 0.2];
    non_centrality_parameter_vec = []; var_vec = []; 
    [~, two_ind] = min(abs(grr_vec-2));
    test_type = 'case-only'; % 'single-locus'; % 'case-controls'; %  % 'case-controls'; % 'case-only';
    test_stat =  'LRT_analytic';  % 'chi-square-analytic'; % % 'chi-square-analytic'; %analytic'; % 'LRT_analytic';
    for i=1:length(mu_vec)
        for j=1:length(f_vec)
            run_j = j 
            for k=1:length(grr_vec)
                p_mat = genetic_relative_risk_to_p_z_x_marginal(f_vec(j), grr_vec(k), mu_vec(i));
                [power_mat, p_vals_vec, stat_vec non_centrality_parameter_vec{i,j}(k)] = ...
                    compute_association_power(p_mat, 1, 1, 0.001, [], test_type, test_stat, [], [], []);
                non_centrality_parameter_vec{i,j}(k) = non_centrality_parameter_vec{i,j}(k)+1; 
                var_vec{i,j}(k) = genetic_relative_risk_to_variance_explained(f_vec(j), grr_vec(k), mu_vec(i), 'binary', 'binary');

            end % loop on grr
        end % loop on allele-frequency
    end % loop on prevalence
    
    figure; hold on; ctr=1;
    for i=1:length(mu_vec)
        for j=1:length(f_vec)
            plot(var_vec{i,j}, non_centrality_parameter_vec{i,j}, color_vec(ctr));
            legend_vec{ctr} = ['\mu=' num2str(mu_vec(i)) '; f=' num2str(f_vec(j))];
            ctr=ctr+1;
        end        
    end
    ctr=1; 
    for i=1:length(mu_vec)
        for j=1:length(f_vec)
            plot(var_vec{i,j}(two_ind), non_centrality_parameter_vec{i,j}(two_ind), [color_vec(ctr) '*']);
            ctr=ctr+1;
        end        
    end
    
    xlabel('Var. explained'); ylabel('E [ LODs]');
    legend(legend_vec, 2); legend('boxoff'); xlim([0 0.05]); 
    my_saveas(gcf, fullfile(new_figs_dir, 'power', 'var_explained_vs_power'), {'epsc', 'pdf'}); 
    
    
    
    
end


total_test_time = cputime - total_test_time






%
% % Small test of how grr is dilluted
% grr = 5; fff = 0.01; mumu = 0.05; alpha=0.1;
%
% p_x_z = genetic_relative_risk_to_p_z_x_marginal(fff, grr, mumu)
% p_x_z_0 = genetic_relative_risk_to_p_z_x_marginal(fff, 1, mumu)
%
% p_x_z_dilluted = alpha .* p_x_z + (1-alpha) .* p_x_z_0;
% [~, grr_dilluted, ~]  = p_z_x_marginal_to_genetic_relative_risk(p_x_z_dilluted)
%
% (grr_dilluted-1) / (grr-1)


