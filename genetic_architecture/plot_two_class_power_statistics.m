% Plot power and N_{1/2} for two-class model
%
% Input:
% power_struct - structure with power for each simulation
% num_samples_struct -  structure with number of samples for each simulation
% f_rare -
% f_rare_vec - vector of rare alleles frequencies
% f_rare_power_vec -
% beta_vec - vector of effect sizes
% N_half_mat - matrix of sample size required to achieve 50% power.
% N_half_mat_only_null_enrichment -
% s_null_vec -  vector of selection coefficients
% show_s_null_ind - some plots only show specific values of selection coefficient s
% i_beta - index of effect size
% i_n - index of ...
% title_str - name of ...
% figs_dir - where to save figures
% new_figs_dir - where to save more figures
% n_samples_half_power_vec -
% power_mat -  matrix of power ...
% n_vec -
% new_kept_alleles_vec -
% alpha_vec -  ??
% new_alpha_vec - ??
% rare_cumulative_per_gene - total rate of mutations in one gene
%  demographic_models_struct - new! allow different demographic models (not just equilibrium)
%
% Output: None. We just plot figures
%
function plot_two_class_power_statistics(equilibrium_parameters_output_file, ...
    power_output_dir, new_figs_dir)



% % % % % num_samples_struct, power_struct, ...
% % % % %     two_class_stat_struct, w_x_null_mat, w_x_harmless, w_all, ...
% % % % %     frac_null_by_freq_cumulative, ...
% % % % %     f_rare, f_rare_vec, beta_vec, ...
% % % % %     s_null_vec, show_s_null_ind,  i_beta, i_n, ...
% % % % %     title_str, figs_dir, new_figs_dir, ...
% % % % %     n_vec, new_kept_alleles_vec, alpha_vec, new_alpha_vec, ...
% % % % %     N, prevalence, gene_length, mu, ...
% % % % %     rare_cumulative_per_gene,  demographic_models_struct) % plot power curves


%%%%% OLDD:
%%%%%%%%%%%%%%
% % % % f_rare, f_rare_vec, f_rare_power_vec, beta_vec, ...
% % % %     N_half_mat, N_half_mat_only_null_enrichment, ...
% % % %     s_null_vec, show_s_null_ind,  i_beta, i_n, ...
% % % %     title_str, figs_dir, new_figs_dir, n_samples_half_power_vec, power_mat,  power_struct{3}.power_mat, ...
% % % %     num_samples_struct{2}.num_samples_mat, ... % num_samples_struct{2}.num_samples_mat2, ...
% % % %     n_vec, new_kept_alleles_vec,  alpha_vec, new_alpha_vec, ...
% % % %     rare_cumulative_per_gene)


AssignGeneralConstants;

orange = [1 0.6 0.1];
population_color_vec = {'k', orange, 'r', 'g', 'b', 'c'};
eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
my_symbol_vec = {'--', '-'}; % flip ordering (set integer powers as solid lines)
selection_color_vec = {'k', 'b', 'g', orange, 'r'}; % replace yellow with orange


load(equilibrium_parameters_output_file);
demographic_models_struct.model_str{end+1} = 'equilibrium';
demographic_models_struct.num_models = length(demographic_models_struct.model_str);


f_vec = 0.001; grr_vec = beta_to_genetic_relative_risk(beta_vec, f_vec, prevalence);
[~, i_beta] = min(abs(grr_vec  - 5)); % 61; % choose specific beta. We choose grr to be 5
[~, i_n] = min(abs(n_vec - 50000)); %  25; % choose specific n
[~, f_ind] = min(abs(f_rare - f_rare_vec))
s_ind = show_s_null_ind(1);


% New: add a power text file output
show_s_null = s_null_vec(show_s_null_ind);
cutoff_freq_vec = [0.01 0.001];

R_all_models = cell(9,1);

n_null_vec_equilibrium = 2 .* N .* mu .* gene_length .* 2 .* ...
    two_class_stat_struct{1}.normalization_factor_x; % Get equilibrium analytic values


for i_d = 2:2 % just the new data 1:demographic_models_struct.num_models
    try
        load(fullfile( power_output_dir, demographic_models_struct.model_str{i_d}, 'two_class_power_parameters.mat')); % load power struct specific for this model
    catch
        continue;
    end
    %    load(power_parameters_output_file);
    
    R_power = cell(20,  length(show_s_null_ind)+1);
    R_power{1,1} = ['Model: ' demographic_models_struct.model_str{i_d}];
    R_power{1,2} = ['\mu_G=' num2str(mu*gene_length)];
    R_power{1,3} = ['\alpha_{birth}=' num2str(alpha_vec)];
    
    ctr=4;     samp_ind = 6;
    R_power{ctr,1} = 'Proportion \rho(T)'; ctr=ctr+1;
    switch demographic_models_struct.model_str{i_d}
        case 'equilibrium'
            rho_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:); % fraction of null alleles conditional on allele freq. below f
            f_null_vec = c_cumulative{1}(show_s_null_ind,:); % fraction of null alleles which are below freq. f
            use_f_vec = f_rare_vec;
            demographic_models_struct.data{i_d}.x_vec = f_rare_vec;
            n_null_vec = 2 .* N .* mu .* gene_length .* 2 .* ...
                two_class_stat_struct{1}.normalization_factor_x; % Get equilibrium analytic values
        otherwise % Schaffner's simulations
            demographic_models_struct.data{i_d} = load(demographic_models_struct.file_names{i_d}); % this gives both s-vec and cumulatives
            rho_vec = demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) ./ ...
                ( demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) + ...
                repmat(demographic_models_struct.data{i_d}.p_vec(end,:), 10, 1) ); % This should go DOWN!!!
            f_null_vec = bsxfun(@rdivide, demographic_models_struct.data{i_d}.p_vec(end:-1:1,:), ...
                demographic_models_struct.data{i_d}.p_vec(end:-1:1,end));  % make sure this is normalized. This should go UP!
            use_f_vec = demographic_models_struct.data{i_d}.x_vec;
            n_null_vec = demographic_models_struct.data{i_d}.p_vec(end:-1:1,end); % take cumulatives
            
    end
    cutoff_freq_inds = zeros(length(cutoff_freq_vec), 1);
    for i_c = 1:length(cutoff_freq_vec)  % loop on possible frequency cutoffs
        [~, cutoff_freq_inds(i_c)] = min( abs(use_f_vec - cutoff_freq_vec(i_c)) );
        R_power{ctr,1} = ['T=' num2str(cutoff_freq_vec(i_c)*100) '%'];
        
        for j=1:length(show_s_null) % loop on selection coefficients
            if(i_c==1)
                R_power{ctr-2,j+1} = ['s=10^(' num2str(round(2*log10(show_s_null(j)))/2) ')'];  % Get headers
            end
            R_power{ctr,j+1} = num2str(rho_vec(j,cutoff_freq_inds(i_c)), 3); % get rho: prob. (null | <= f)
        end
        ctr=ctr+1;
    end
    R_power{ctr,1} = 'T_{optimal}';
    for j=1:length(show_s_null) % Add T_opt
        [T_min_num_samples(j) T_min_ind(j)] = min( round(0.5*num_samples_struct{samp_ind}.num_samples_mat(1:end-0, j)) ); % show_s_null_ind(j))); % get threshold achieving maximal power (minimal sample size)
        [T_min_num_samples_eric(j) T_min_ind_eric(j)] = min( num_samples_struct{samp_ind+3}.num_samples_mat(1:end-0, j) ); % show_s_null_ind(j))); % get threshold achieving maximal power (minimal sample size)
        R_power{ctr,j+1} = num2str(rho_vec(j, T_min_ind(j)), 3);
    end
    ctr=ctr+2;
    R_power{ctr,1} = 'Total Allele Freq. n_{null}(1)'; ctr=ctr+1;
    for j=1:length(show_s_null) % loop on selection coefficients
        R_power{ctr,j+1} = num2str(n_null_vec(j), 3);
    end
    
    ctr=ctr+2;
    R_power{ctr,1} = 'Cumulative Allele Dist. f_{null}(T)'; ctr=ctr+1;
    for i_c = 1:length(cutoff_freq_vec)  % loop on possible frequency cutoffs
        R_power{ctr,1} = ['T=' num2str(cutoff_freq_vec(i_c)*100) '%'];
        for j=1:length(show_s_null) % loop on selection coefficients
            R_power{ctr,j+1} = num2str(f_null_vec(j, cutoff_freq_inds(i_c)), 3); % ['s=10^(' num2str(log10(s_null_vec(j))) ')'];
        end
        ctr=ctr+1;
    end
    R_power{ctr,1} = 'T_{optimal}';
    for j=1:length(show_s_null) % Add T_opt
        R_power{ctr,j+1} = num2str(f_null_vec(j, T_min_ind(j)), 3);
    end
    
    
    ctr=ctr+2;
    R_power{ctr,1} = 'Apparent enrichment for true GRR=5'; ctr=ctr+1;
    for i_c = 1:length(cutoff_freq_vec)  % loop on possible frequency cutoffs
        R_power{ctr,1} = ['T=' num2str(cutoff_freq_vec(i_c)*100) '%'];
        for j=1:length(show_s_null) % loop on selection coefficients
            R_power{ctr,j+1} = num2str(1 + (5-1) * rho_vec(j,cutoff_freq_inds(i_c)), 3); % ['s=10^(' num2str(log10(s_null_vec(j))) ')'];
        end
        ctr=ctr+1;
    end
    R_power{ctr,1} = 'T_{optimal}';
    for j=1:length(show_s_null) % Add T_opt
        [T_min_num_samples T_min_ind] = min( round(0.5*num_samples_struct{samp_ind}.num_samples_mat(1:end-0, j)) ); % show_s_null_ind(j))); % get threshold achieving maximal power (minimal sample size)
        [T_min_num_samples_eric T_min_ind_eric] = min( num_samples_struct{samp_ind+3}.num_samples_mat(1:end-0, j) ); % show_s_null_ind(j))); % get threshold achieving maximal power (minimal sample size)
        R_power{ctr,j+1} = num2str(1 + (5-1) * rho_vec(j,T_min_ind), 3);
    end
    
    
    ctr=ctr+2;
    R_power{ctr,1} = 'Required Sample Size (50% power, pval<2.5*10^6, GRR=5, prev.=5%). Num. cases required assuming either:  Chi-square, equal #controls, or: (Eric''s Zscore, with inifinite #controls)'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (missense, perfect knowledge)'; ctr=ctr+1;  % add missense and LOF separately
    R_power{ctr,1} = 'SS (missense<= 1%)'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (missense<= 0.1%)'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (missense<= T_{optimal})'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (all LOF, perfect knowledge)'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (all LOF + missense<= 1%)'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (all LOF + missense<= 0.1%)'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (all LOF + missense<= T_{optimal})'; ctr=ctr+1;
    R_power{ctr,1} = 'T_{optimal}'; ctr=ctr+1;
    R_power{ctr,1} = 'SS (ratio:  missense<= T_{optimal} / perfect knowledge)'; ctr=ctr-9;
    
    for j=1:length(show_s_null) % loop on selection coefficients and compute power
        j_is = j
        i_d_is = i_d
        % num_samples_struct{6} % here we have mixture of nulls and non-nulls
        % num_samples_struct{8} % here we have pure nulls
        
        R_power{ctr,j+1} = [num2str(round(0.5*num_samples_struct{samp_ind+1}.num_samples_mat(end, j))) ...  % sample size for pure missense
            ' (' num2str(num_samples_struct{samp_ind+4}.num_samples_mat(end, j)) ')'];  % show_s_null_ind(j)))
        R_power{ctr+4,j+1} = [num2str(round(0.5*num_samples_struct{samp_ind+2}.num_samples_mat(end, j))) ...  % sample size for pure missense
            ' (' num2str(num_samples_struct{samp_ind+5}.num_samples_mat(end, j)) ')']; ctr=ctr+1;
        for i_c = 1:length(cutoff_freq_vec)  % loop on possible frequency cutoffs
            i_c_is = i_c
            [~, cur_f_ind] = min( abs( use_f_vec - cutoff_freq_vec(i_c) ) ); % find index in
            R_power{ctr,j+1} = [num2str(round(0.5*num_samples_struct{samp_ind}.num_samples_mat(cur_f_ind, j))) ... % sample size for dilluted missense, different cutoffs
                ' (' num2str(num_samples_struct{samp_ind+3}.num_samples_mat(cur_f_ind, j)) ')']; % show_s_null_ind(j)));
            R_power{ctr+4,j+1} = [num2str(round(harmmean([0.5*num_samples_struct{samp_ind}.num_samples_mat(cur_f_ind, j) ...
                0.5*num_samples_struct{samp_ind+2}.num_samples_mat(end, j)])/2)) ...
                ' (' num2str(round(harmmean([ 0.5*num_samples_struct{samp_ind+3}.num_samples_mat(cur_f_ind, j) ...
                num_samples_struct{samp_ind+5}.num_samples_mat(end, j)])/2)) ')'];  ctr=ctr+1;
        end
        [T_min_num_samples T_min_ind] = min( round(0.5*num_samples_struct{samp_ind}.num_samples_mat(1:end-0, j)) ); % show_s_null_ind(j))); % get threshold achieving maximal power (minimal sample size)
        [T_min_num_samples_eric T_min_ind_eric] = min( num_samples_struct{samp_ind+3}.num_samples_mat(1:end-0, j) ); % show_s_null_ind(j))); % get threshold achieving maximal power (minimal sample size)
        R_power{ctr,j+1} = [num2str(T_min_num_samples) ' (' num2str(T_min_num_samples_eric)  ')']; % sample size for dilluted, optimal threshold
        R_power{ctr+4,j+1} = [num2str( round(harmmean([T_min_num_samples 0.5*num_samples_struct{samp_ind+2}.num_samples_mat(end, j)])/2) ) ...
            ' (' num2str( round(harmmean([T_min_num_samples_eric num_samples_struct{samp_ind+5}.num_samples_mat(end, j)])/2) ) ')'];
        ctr=ctr+1; % num samples when taking optimal threshold
        R_power{ctr+4,j+1} = [num2str(100*use_f_vec(T_min_ind), 3) '% (' ...
            num2str(100*use_f_vec(T_min_ind_eric), 3) '%)']; % took optimal threshold. Should be: use_f_vec
        R_power{ctr+5,j+1} = num2str(str2nums(R_power{ctr-1,j+1}) ./ str2nums(R_power{ctr-4,j+1}), 3); ctr=ctr-4;
        
        
        %         for combined_ctr = 1:4; % here we combine tests. Sample size (approximate) is harmonic mean
        %             R_power{ctr+4+combined_ctr,j} = harmmean(R_power{ctr+4,j+1}, R_power{ctr+combined_ctr,j+1})/2;
        %
        %         end
        
    end
    R_power = num2str_cell(R_power, 3);
    savecellfile(num2str_cell(R_power), fullfile(new_figs_dir, demographic_models_struct.model_str{i_d}, 'power_table.txt'), ...
        [], 1);
    
    power_str_ind = 22;
    for i_plot = 1:9 % Copy all to unique structure for all models
        R_all_models{i_plot}{1,1} = R_power{power_str_ind+i_plot,1};
        R_all_models{i_plot}{2,1} = 'Model:';
        R_all_models{i_plot}{i_d+2,1} = demographic_models_struct.model_str{i_d};
        for j=1:length(show_s_null)
            R_all_models{i_plot}{2,j+1} = ['s=10^(' num2str(round(2*log10(show_s_null(j)))/2) ')'];  % Get headers
            R_all_models{i_plot}{i_d+2, j+1} = R_power{power_str_ind+i_plot, j+1};
        end
    end
    
    
    
    %%%%%%%%%%% End preparing table %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % New: Plot contour plots for equal sample size as a function of selection coefficient and effect size
    
    % for i_d = 1:demographic_models_struct.num_models
    stat_ctr=1;
    for stat_used = {'chi-square-analytic', 'eric_crude_enrichment_analytic'}
        
        figure; [C,h] = contour(s_null_vec(show_s_null_ind), grr_vec, ...
            log10(num_samples_struct{4+stat_ctr-1}.num_samples_mat(1:10,:)'), 1:0.5:10, ...
            'linewidth', 2); % 5 is for my power statistic, 6 is for Eric's crude statistic
        % [C,h] = contour(new_kept_alleles_vec, new_alpha_vec, power_struct{3}.power_mat);
        set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
        th = clabel(C,h) % change text to include powers of 10
        for i=1:0.5:10
            set(findobj(th,'String',num2str(i)),'String',['10^{' num2str(i) '}'])
        end
        set(gca,'xscale','log'); % set(gca, 'yscale', 'log');
        colormap cool;
        xlabel('s'); ylabel('GRR'); title(['Isocurves for sample size required for 90% power. \alpha=' ...
            num2str(num_samples_struct{4}.p_val_cutoff) '. \pi=' num2str(prevalence)]);
        ylim([1 10]);
        add_faint_grid(0.5);
        
        % % % figure; my_contour(log10(s_null_vec), grr_vec, log10(num_samples_struct{4}.num_samples_mat'), 1:0.5:5);
        % % % ylim([1 10]); legend( num2str(  (1:0.5:5)' ) )
        
%         my_saveas(gcf, fullfile(new_figs_dir, demographic_models_struct.model_str{i_d}, ...
%             ['fig1c_sample_size_gene_burden_test_function_of_GRR_and_s_isocurves_' stat_used{1}]), ...
%             {'epsc', 'pdf'});
        my_saveas(gcf, fullfile(new_figs_dir, 'power/isocurves', ...
            [demographic_models_struct.model_str{i_d} '_sample_size_isocurves_' stat_used{1}]), ...
            {'epsc', 'pdf'});
        
        

        
        % NEW! Plot just power for different values of GRR
        figure; show_grr_vec = [1.3 1.5 2 3 5 10];
        if( strcmp(demographic_models_struct.model_str{i_d}, 'equil') && strcmp(stat_used{1}, 'eric_crude_enrichment_analytic') ) % plot also power by theory
            for j=1:length(show_grr_vec)
                
                tmp_eric_sample_size = 2.*( norminv(p_val_cutoff) - sqrt(show_grr_vec(j)) .* norminv(0.9) ).^2 ./ ...
                    ((show_grr_vec(j)-1).^2 .* 0.4 .* vec2row(two_class_stat_struct{1}.normalization_factor_x));
                for i_s = 1:length(s_null_vec)
                    p_mat = zeros(1,4);
                    p_mat(4) = show_grr_vec(j) * prevalence *  0.4 .* two_class_stat_struct{1}.normalization_factor_x(i_s);
                    p_mat(2) = prevalence - p_mat(4);
                    p_mat(3) =  0.4 .* two_class_stat_struct{1}.normalization_factor_x(i_s) - p_mat(4);
                    p_mat(1) = 1-sum(p_mat);
                    if(min(p_mat) < 0)
                        tmp_eric_sample_size(i_s)  = -1;
                    end
                    %                    p_mat = genetic_relative_risk_to_p_z_x_marginal(0.4 .* two_class_stat_struct{1}.normalization_factor_x(i_s), show_grr_vec(j), prevalence);
                    %                    show_rr_vec(j,i_s) = p_mat(:,4) ./ ( sum(p_mat(:, [2 4])) .* sum(p_mat(:, [3 4])) );
                end
                loglog(s_null_vec, tmp_eric_sample_size, color_vec(j), 'linewidth', 2); hold on; % [color_vec(j) '--']
            end
        else % here a traditional curve
            for j=1:length(show_grr_vec)
                [~, grr_ind] = min(abs(grr_vec - show_grr_vec(j)))
                loglog(show_s_null, num_samples_struct{4+stat_ctr-1}.num_samples_mat(1:10,grr_ind), ...  % no need to flip s
                    color_vec(j), 'linewidth', 2); hold on;
            end
        end % if strcmp 

        
        xlabel('s'); ylabel('Sample Size'); 
        % title(['Sample size required for 90% power. \alpha=' ...
        %    num2str(num_samples_struct{4}.p_val_cutoff) '. \pi=' num2str(prevalence) '. ' ...
        %    demographic_models_struct.model_str{i_d}]);
        add_faint_grid(0.5);
        h_leg = legend(cellstr([repmat('RR=', 6, 1) num2str(show_grr_vec')]), 2);
        set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
        xlim([10^(-5) 0.1]); ylim([100 10^7]); 

% % Show n_null on top axis: 
%         ax2=axes('Position',get(gca,'Position'),...
%             'XAxisLocation','top',...
%             'YAxisLocation','right',...
%             'Color','none',...
%             'XColor','k','YColor',[0.8 0.8 0.8]);
%         x_labels = n_null_vec_equilibrium(show_s_null_ind([1 2:2:end])); x_labels(1) = 0.4;
%         line(-n_null_vec_equilibrium, zeros(length(n_null_vec_equilibrium), 1),'Color',[0.8 0.8 0.8],'Parent',ax2);
%         set(ax2,'xscale','log');
%         set(ax2, 'xtick', [-1 -x_labels(2:end-1)' -0.00012]); 
%         set(ax2, 'xticklabel', cellstr(num2str(x_labels, 2)));
%         xlabel('n_{null}');        
        
        my_saveas(gcf, fullfile(new_figs_dir, 'power/', ...
            [demographic_models_struct.model_str{i_d} '_sample_size_' stat_used{1}]), ...
            {'epsc', 'pdf'});
        
        
        stat_ctr=stat_ctr+1;
    end % loop on stat counter
    
    
    demographic_models_struct.data{i_d}.sample_size_inflation_factor =  1 ./ (rho_vec .* f_null_vec);
    
    continue; % for now don't do other plots below
    
    
    
    ctr=5;
    
    figure; % Plot power for different values of s
    y_lim = [min( min(min(num_samples_struct{ctr+1}.num_samples_mat(:,:), ... % show_s_null_ind))), ...
        min(min(num_samples_struct{ctr+2}.num_samples_mat(:,:))) ), ... % show_s_null_ind))) ) ...
        max( max(max(num_samples_struct{ctr+1}.num_samples_mat(:,:))), ... % show_s_null_ind))), ...
        max(max(num_samples_struct{ctr+2}.num_samples_mat(:,:))) ) ) )]; % show_s_null_ind))) )];
    y_lim(1) = 10^(floor(log10(y_lim(1))));
    y_lim(2) = 10^(ceil(log10(y_lim(2))));
    for i_p=1:2 % plot N_1/2 for different values of s and f^*
        for log_x_flag = 1:1 % 0:1 % 1:1 % plot only on log scale
            subplot(1,2, i_p);
            %            subplot(2,2, i_p + log_x_flag*2);
            switch log_x_flag
                case 0
                    log_x_str = '';
                    semilogy(f_rare_vec, beta_vec(i_beta)^2 .* 0.5 .*(num_samples_struct{ctr+i_p,q}.num_samples_mat(:,:)), ... % show_s_null_ind)), ...
                        'linewidth', 2); hold on; % 0.5 is because we look at people, not chromosomes
                    semilogy(f_rare_vec, beta_vec(i_beta)^2 .* 0.5 .*(num_samples_struct{ctr+2+i_p}.num_samples_mat(:,:)), ... % show_s_null_ind)), ...
                        '--', 'linewidth', 2); % 0.5 is because we look at people, not chromosomes
                    %                 semilogy(f_rare_vec, beta_vec(i_beta)^2 .* 0.5 .*(N_half_mat_only_null_enrichment2{i_p}(:,show_s_null_ind)), ...
                    %                     ':', 'linewidth', 2); % 0.5 is because we look at people, not chromosomes
                    
                case 1
                    log_x_str = ' (log)';
                    loglog(f_rare_vec, beta_vec(i_beta)^2 .* 0.5 .*(num_samples_struct{ctr+i_p}.num_samples_mat(:,:)), ... % show_s_null_ind)), ...
                        'linewidth', 2); hold on; % 0.5 is because we look at people, not chromosomes
                    loglog(f_rare_vec, beta_vec(i_beta)^2 .* 0.5 .*(num_samples_struct{ctr+2+i_p}.num_samples_mat(:,:)), ... % show_s_null_ind)), ...
                        '--', 'linewidth', 2); % 0.5 is because we look at people, not chromosomes
                    %                 loglog(f_rare_vec, beta_vec(i_beta)^2 .* 0.5 .*(N_half_mat_only_null_enrichment2{i_p}(:,show_s_null_ind)), ...
                    %                     ':', 'linewidth', 2); % 0.5 is because we look at people, not chromosomes
            end
            xlabel(['allele freq. cutoff' log_x_str]); ylabel('Sample size (1/2 power) X \beta^2 (log)');
            xlim([0 0.95]);
            cur_title_str = strrep(title_str, 'Detect. power rare.', 'N_{1/2} \times \beta^2, ')
            cur_title_str = strrep(cur_title_str, 'f^*=1%,', '');
            cur_title_str = strdiff(cur_title_str, 'cum. freq.,');
            %            title(cur_title_str);
        end % loop on log flag
        
        % %         figure; % Make contour plot of N_{1/2}
        % %         [C,h] = contour(s_null_vec, f_rare_vec, N_half_mat);
        % %         set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
        % %         colormap cool;
        % %         xlabel('Selection s'); ylabel('Allele freq. f^*');
        % %         title('Sample size needed to acheieve 1/2 power');
        % %         my_saveas(gcf, fullfile(figs_dir, 'power', 'rare_variants_power_contour_as_function_of_effect_size'), ...
        % %         {'epsc', 'pdf'});
        %    y_lim = get(gca, 'ylim');
        if(y_lim(2) > y_lim(2))
            ylim(y_lim); % temp !!!
        end
    end % loop on p-value
    cur_title_str = ['N_{1/2} \times \beta^2,  \alpha=33.3% null, ' ...
        'c=100%, p-val cutoff =  0.5 (left), 2.5e-006 (right)'];
    %    suptitle(cur_title_str);
    legend_vec = [repmat('s^*= -', length(show_s_null_ind), 1) num2str(vec2column(s_null_vec(show_s_null_ind)), 2)];
    
    
    legend(legend_vec,1); % 3);
    my_saveas(gcf, fullfile(new_figs_dir, 'power', ...
        'rare_variants_sample_size_N_0.5_as_function_of_s_and_f_star'), {'epsc', 'pdf'});
    
    %%%%%%%%%Plot Ratio of sample size needed (enrichment vs. no enrichment) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    ratio_vec =  num_samples_struct{ctr+2}.num_samples_mat(:,show_s_null_ind) ./  num_samples_struct{ctr+2+2}.num_samples_mat(:,show_s_null_ind);
    [min_sample_size min_sample_size_inds] = min(num_samples_struct{ctr+2}.num_samples_mat(:,show_s_null_ind));
    loglog(f_rare_vec, ratio_vec, 'linewidth', 2); hold on; % 0.5 is because we look at people, not chromosomes
    plot(f_rare_vec(min_sample_size_inds), diag(ratio_vec(min_sample_size_inds,:)), '*k', 'markersize', 8);
    xlabel(['allele freq. cutoff' log_x_str]); ylabel('Fold-Saving by Enrichment');
    legend(legend_vec,2); % 3);
    my_saveas(gcf, fullfile(new_figs_dir, 'power', ...
        'rare_variants_sample_size_N_0.5_enrichment_effect_as_function_of_s_and_f_star'), {'epsc', 'pdf'});
    
    optimal_ratio_vec = min(num_samples_struct{ctr+2}.num_samples_mat(:,show_s_null_ind)) ./ min(num_samples_struct{ctr+2+2}.num_samples_mat(:,show_s_null_ind));
    figure; semilogx(s_null_vec(show_s_null_ind), optimal_ratio_vec, '*');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot sample size as function of effect size beta 
    figure; plot(beta_vec, num_samples_struct{1}.num_samples_mat); title('# samples needed as function of \beta');
    xlabel('\beta'); ylabel('N_{1/2}');
    
    
    beta_vec_str = beta_vec-1; % string to display
    figure; hold on; imagesc_with_labels(power_struct{1}.power_mat(end:-1:1,:), n_vec, beta_vec(end:-1:1), 0.5, [], 10);     % Plot power
    ylabel('Effect size \beta'); xlabel('Sample size');
    title(title_str);
    
    figure; % Make contour plot of power. Plot equi-power curves as function of samples and effect size plain
    [C,h] = contour((n_vec(1:51)), (beta_vec+1), power_struct{1}.power_mat(:,1:51));
    set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
    colormap cool;
    xlabel('Sample size (n)'); ylabel('Effect size (\beta)');
    title(title_str);
    my_saveas(gcf, fullfile(figs_dir, 'power', 'rare_variants_power_contour_as_function_of_effect_size'), ...
        {'epsc', 'pdf'});
    
    % Plot power as function of CAF f  
    figure; plot(f_rare_vec, power_struct{2}.power_mat, 'linewidth', 2);
    title([ strdiff(title_str,  ['f^*=' num2str(f_rare*100,2) '%, ']) ', n=' num2str(n_vec(i_n))]);
    xlabel('allele freq. cutoff f^*'); ylabel('Power');
    my_saveas(gcf, fullfile(figs_dir, 'power', 'rare_variants_power_as_function_of_allele_freq_cutoff'), ...
        {'epsc', 'pdf'});
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute and plot power when enriching for functional 'null' alleles.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure; hold on; % Make contour plot of power
    [C,h] = contour(new_kept_alleles_vec, new_alpha_vec, power_struct{3}.power_mat);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
    colormap cool;
    plot(new_kept_alleles_vec, alpha_vec ./ new_kept_alleles_vec, '--', 'linewidth', 3); % plot feasability boundary
    plot(1, alpha_vec, '*', 'markersize', 20); % plot current scenario
    plot(alpha_vec, 1, '*r', 'markersize', 20); % plot best scenario
    xlabel('Fraction of alleles retained (r)'); ylabel('Fraction of null alleles out of retained');
    
    title(['Power as function of enrichment \alpha and fraction of retained alleles r. Sample size=' ...
        num2str(n_vec(i_n)) ', c=' num2str(rare_cumulative_per_gene*100,3) '%, \beta=' num2str(beta_vec(i_beta), 3)]);
    legend('power', 'best feasible curve', 'all alleles', 'only functional alleles');
    my_saveas(gcf, fullfile(figs_dir, 'power', 'rare_variants_power_contour_as_function_of_enrichment_alpha'), ...
        {'epsc', 'pdf'});
    
    
    figure; hold on; % Make contour plot of N_{1/2}, as function of ...
    [C,h] = contour(new_kept_alleles_vec, new_alpha_vec, log10(num_samples_struct{2}.num_samples_mat));
    set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
    colormap cool;
    plot(new_kept_alleles_vec, alpha_vec ./ new_kept_alleles_vec, '--', 'linewidth', 3); % plot feasability boundary
    plot(new_kept_alleles_vec, repmat(alpha_vec, length(new_kept_alleles_vec), 1), '--g', 'linewidth', 3); % plot feasability boundary
    plot(1, alpha_vec, '*', 'markersize', 20); % plot current scenario
    plot(alpha_vec, 1, '*r', 'markersize', 20); % plot best scenario
    xlabel('Fraction of alleles retained (r)'); ylabel('Fraction of null alleles out of retained');
    %title(['N_{1/2} as function of enrichment \alpha and fraction of retained alleles r. Sample size=' ...
    %    num2str(n_vec(i_n)) ', c=' num2str(rare_cumulative_per_gene*100,3) '%, \beta=' num2str(beta_vec(i_beta), 3)]);
    legend('N_{1/2}', 'best feasible curve', 'lower-bound curve', 'all alleles', 'only functional alleles');
    my_saveas(gcf, fullfile(figs_dir, 'power', 'rare_variants_N_50_contour_as_function_of_enrichment_alpha'), ...
        {'epsc', 'pdf'});
    
end % loop on demographic models

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure; % New: plot inflation factor as function of threshold f for different values of s and different demographic models
for j=1:length(show_s_null)
    for include_LOF_flag = 0:1 % plot I_T or I_star_T (latter includes also LOFs)
        figure; model_ctr=1; model_legend = [];
        for i_d = 1:demographic_models_struct.num_models
            if(ismember(demographic_models_struct.model_str{i_d}, {'equil', 'expan1', 'expan2', 'europ', '2phase'}))
                if(include_LOF_flag)
                    cur_inflation_factor_vec = ...
                        harmmean([(5/4).*demographic_models_struct.data{i_d}.sample_size_inflation_factor(j,:);  ...
                        repmat(5, 1, length(demographic_models_struct.data{i_d}.x_vec))])./2;
                    inflation_str = '_star'; inflation_str_numeric = '^*';
                else
                    cur_inflation_factor_vec = (5/4).*demographic_models_struct.data{i_d}.sample_size_inflation_factor(j,:);
                    inflation_str = ''; inflation_str_numeric = '';
                end
                
                semilogx(demographic_models_struct.data{i_d}.x_vec, cur_inflation_factor_vec, ...
                    color_vec(model_ctr), 'linewidth', 2);  hold on;
                model_legend = [model_legend {demographic_models_struct.model_str{i_d}}];
                model_ctr = model_ctr+1;
            end
        end
        if(~include_LOF_flag)
            ylim([1 10]);
        else
            ylim([1 5.01]);
        end
        xlabel('Derived Allele Freq.'); ylabel('Sample Size Inflation Factor');
        title(['Inflation $I' inflation_str_numeric '(T)$ in Sample size as function of thereshold. $s=10^{' ...
            num2str(round(2*log10(show_s_null(j)))/2) '}$'], 'interpreter', 'latex');
        %     legend_vec = [repmat('s= 10^{', length(show_s_null), 1) num2str(log10(show_s_null'),3) ...
        %         repmat('}', length(show_s_null), 1)];
        %     legend_vec = cellstr(legend_vec); legend_vec = strrep_cell(legend_vec, ' ', ''); legend_vec{1} = 's= 0'; % fix s=0
        legend(model_legend); legend ('boxoff');
        my_saveas(gcf, fullfile(new_figs_dir, 'all_models', ['I' inflation_str '_T'], ...
            strrep(['power_inflation_factor_I' inflation_str '_T_s_10_'  num2str(round(2*log10(show_s_null(j)))/2)], '.', '_' )), ...
            {'epsc', 'pdf'});
    end % loop on include LOF flag
end % loop on selection coefficient


for i_plot=1:9
    save_i = i_plot
    cur_save_file = fullfile(new_figs_dir, 'all_models', ...
        [strrep(strrep(strrep( strrep(strrep(strrep(strrep(R_all_models{i_plot}{1,1}, ' ', '_'), '(', '_'), ')', '_'), '%', ' pct'), '<=', '_less'), '/', '_'), ':', '_'),  ...
        '_power_table.txt'])
    savecellfile(num2str_cell(R_all_models{i_plot}), cur_save_file,  [], 1);
    
    if(i_plot==1) % Also plot 'pure' sample size
        cur_plot_mat = cell2mat(str2nums_cell(R_all_models{i_plot}(3:end,2:end), 1));
        
        cur_model_legend = R_all_models{i_plot}(2+find(~isempty_cell(R_all_models{i_plot}(3:end,1))),1);

        bad_inds = strmatch('varsel', cur_model_legend); good_inds = setdiff(1:length(cur_model_legend), bad_inds); % remove variable selection 
        cur_model_legend = cur_model_legend(good_inds); 
        cur_plot_mat = cur_plot_mat(good_inds,:); 
        
        
        %     for i=1:size(cur_plot_mat, 1)
        %         for j=1:size(cur_plot_mat, 2)
        %             cur_plot_mat{i,j} = cur_plot_mat{i,j}(1);
        %         end
        %     end
        %     cur_plot_mat = cell2mat(cur_plot_mat);
        %
        figure;
        for j=1:size(cur_plot_mat, 1)
            switch cur_model_legend{j}
                case {'equil', 'expan1', 'expan2', '2phase'}
                    cur_symbol_vec = '-';
                case {'ice', 'finn1', 'finn2', 'europ'}
                    cur_symbol_vec = '--';
                case {'varsel1', 'varsel2'}
                    cur_symbol_vec = ':';
            end
            loglog(show_s_null, cur_plot_mat(j,1:end)', [color_vec(j) cur_symbol_vec], 'linewidth', 2); hold on; % no need to flip power 
        end
        xlabel('s'); ylabel('Sample size'); %title('Sample Size for different models. Prevalence=0.05. GRR=5. Power=50%');
        h_leg = legend(cur_model_legend, 2); % legend('boxoff');
        set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
        add_faint_grid(0.5);        
        x_ticks = get(gca, 'xtick'); x_tick_labels = get(gca, 'xticklabel');

        ax1=gca;
        %        set(gca, 'xtick', show_s_null); 
        ax2=axes('Position',get(gca,'Position'),...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none',...
            'XColor','k','YColor',[0.8 0.8 0.8]);
        x_labels = n_null_vec_equilibrium(show_s_null_ind([1 2:2:end])); x_labels(1) = 0.4;
        hl2=line(-n_null_vec_equilibrium, repmat(0, length(n_null_vec_equilibrium), 1),'Color',[0.8 0.8 0.8],'Parent',ax2);
        set(ax2,'xscale','log');
        set(ax2, 'xtick', [-1 -x_labels(2:end-1)' -0.00012]); 
        set(ax2, 'xticklabel', cellstr(num2str(x_labels, 2)));
        xlabel('n_{null}');
        
        
        % Add additional x-axis for n_null 
        my_saveas(gcf, fullfile(new_figs_dir, 'all_models', 'sample_size_different_models'), {'epsc', 'pdf'});
    end
    
    
end


selection_legend_vec = [repmat('s= 10^{', length(show_s_null), 1) num2str(log10(show_s_null'),3) ...
    repmat('}', length(show_s_null), 1)];
selection_legend_vec = selection_legend_vec(6:10,:);
selection_legend_vec = strrep_cell(cellstr(selection_legend_vec), ' ', ''); 

figure;
%for i=1:length(power_struct{4}.lambda_vec) % NEW! Plot power as function of sample size - for NHGRI
for s_ind = 6:length(show_s_null)
    semilogx(power_struct{4}.n_vec, squeeze(power_struct{4}.power_mat(10,s_ind,:)), ...
        'color', selection_color_vec{ceil(s_ind/2)}, ...
        'linestyle', my_symbol_vec{mod_max(s_ind,2)}, 'linewidth', 2); hold on;
end
%    xlabel('Number of cases');
ylabel('Power');
h_leg = legend(selection_legend_vec, 2);
set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
xlim([10^(2.0), 3*10^5]); % ylim([-0.2 1]);
%set(gca, 'xtick', [10^2 5*10^2 10^3 5*10^3 10^4 5*10^4 10^5 5*10^5]); 
%set(gca, 'xticklabel', {'', '', '', '', '', '', '', ''});

%{'10^2', '5\times10^2', '10^3', '5\times10^3', '10^4', '5\times10^4', '10^5', '5\times10^5'});  

% Add xlabels
lambda_ctr=1;
text(37, -0.032, ['1+\lambda=' num2str(power_struct{4}.lambda_vec(10))]); % default lambda
%for k=2:5
%   text(0.75*10^k, -0.03, ['10^' num2str(k)]);   
%   text(3.0*10^k, -0.03, ['5\times10^' num2str(k)]);   
%end
lambda_y_plot_vec = [-0.07 -0.11];
for j_lambda = [8 12]
    text(37, lambda_y_plot_vec(lambda_ctr)-0.002, ['1+\lambda=' num2str(power_struct{4}.lambda_vec(j_lambda))]);
    for k=2:5
        tmp_str =  num2str(10^k*g_power(power_struct{4}.lambda_vec(10)-1)/g_power(power_struct{4}.lambda_vec(j_lambda)-1), 2);
        if(isempty(strfind(tmp_str, 'e')))
            if(length(tmp_str) == 2)
               tmp_str = [tmp_str(1) '.' tmp_str(2) '\times10^1']; 
            else
                tmp_str = [tmp_str '\times10^0']; 
            end
        else
        tmp_str = strrep(tmp_str, 'e+', '\times 10^{'); 
        tmp_str = strrep(strrep([tmp_str '}'], '{00', '{'), ' ', ''); 
        end
        tmp_str = strrep(tmp_str, '1\times', '1.0\times'); 
        text(0.75*10^k, lambda_y_plot_vec(lambda_ctr), tmp_str); 
    end
    lambda_ctr=lambda_ctr+1;
end
add_faint_grid(0.5);

my_saveas(gcf, fullfile(new_figs_dir, 'all_models', 'power_vs_sample_size_and_s'), {'epsc', 'pdf'});
