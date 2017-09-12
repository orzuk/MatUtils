% Plot several statistics for the two-class model
%
% Input:
% two_class_stat_struct - structure with parameters
% w_x_null_mat - ??
% w_x_harmless - ??
% w_all - ??
% f_rare_vec - vector of allele frequencies f*
% frac_null_by_freq_cumulative - Matrix of Prob(allele null | freq. <= f). Notation:  \alpha_{s^*}(<= f)
% xi_prob_leq_f - Matrix of probability of allele being below prob. f. Notation: \xi( <= f)
% N - effective population size
% show_s_null - vector of specific values of s
% show_s_null_ind - vector of indices for specific values of s
% s_null_vec - vector of s for the null alleles
% rare_cumulative_per_gene - number of total rare variants per gene (1)
% alpha_vec - proportion of variants which are null
% title_str - string describing the two-class model parameters
% figs_dir - where to save figures
% new_figs_dir - where to save other figures
% figs_for_paper_flag - if this is 'on', just plot the figures needed for paper (if 'off', plot additional figures)
% demographic_models_struct - NEW! WE allow here more population models !!!
%
% Output:
% Compute several statistics for distribution of allele frequencies
%
function plot_two_class_equilibrium_statistics(two_class_stat_struct, w_x_null_mat, w_x_harmless, w_all, ...
    f_rare_vec, frac_null_by_freq_cumulative, xi_prob_leq_f, ...
    N, show_s_null, show_s_null_ind, s_null_vec, rare_cumulative_per_gene, alpha_vec, prevalence, ...
    title_str, figs_dir, new_figs_dir, plot_bayes_factor, figs_for_paper_flag, ...
    demographic_models_struct, equilibrium_parameters_output_file)

AssignGeneralConstants;
eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
my_symbol_vec = {'--', '-'}; % dashed for s=0 and 1/2 powers of 10, solid for whole powers of 10

if(~exist('plot_bayes_factor', 'var') || isempty(plot_bayes_factor))
    plot_bayes_factor = 0; % plot independent of mixture coefficient
end
if(~exist('figs_for_paper_flag', 'var') || isempty(figs_for_paper_flag))
    figs_for_paper_flag = 0;
end
if(~exist('prevalence', 'var') || isempty(prevalence))
    prevalence=0.05; % assume a fixed prevalence (5%)
end
if(~exist('demographic_models_struct', 'var'))
    demographic_models_struct = []; demographic_models_struct.model_str = [];
end
demographic_models_struct.model_str{end+1} = 'equilibrium';
demographic_models_struct.file_names{end+1} = ''; % add one for equilibrium
demographic_models_struct.num_models = length(demographic_models_struct.model_str);

%    legend_vec = [repmat('s^*= -', length(show_s_null), 1) num2str(show_s_null',3)];
legend_vec = [repmat('s= 10^{', length(show_s_null), 1) num2str(log10(show_s_null'),3) ...
    repmat('}', length(show_s_null), 1)];
legend_vec = cellstr(legend_vec);
legend_vec = strrep_cell(legend_vec, ' ', '');
legend_vec{1} = 's= 0'; % fix s=0

prob_null_given_x_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,end);


R_med = cell(demographic_models_struct.num_models+2, length(show_s_null)+2); % New! prepare table
R_med{1,1} = 's';
R_med{3,1} = 'f_{null}(equil.)';

[mu, gene_length, mu_per_gene] = set_genomic_mutation_rates(); % set mutation rates

Big_S = show_s_null .* N .* 4;
normalization_factor_vec = phi_s_integral(0.999999999,-Big_S, 1) - phi_s_integral(0.000000001,-Big_S, 1);

n_null_vec_equilibrium = 2 .* N .* mu .* gene_length .* 2 .* normalization_factor_vec; % Get equilibrium analytic values

for j=1:length(show_s_null)
    R_med{1, j+1} = ['$10^{' num2str(log10(show_s_null(j)),3) '}$'];
    R_med{3, j+1} = num2str(n_null_vec_equilibrium(j), 3);
end
R_med{1,2} = '0';

R_med{5,1} = 'Median Allele Freq.:'; R_ctr=6;
for i_d = 1:demographic_models_struct.num_models
    run_id = i_d
    if(i_d > length(demographic_models_struct.file_names))
        continue;
    end
    if(~ismember(demographic_models_struct.model_str{i_d}, {'equil', '2phase', 'expan1', 'expan2', 'europ', 'ice', 'finn1', 'finn2'}))
        continue;
    end
    demographic_models_struct.data{i_d} = load(demographic_models_struct.file_names{i_d}); % this gives both s-vec and cumulatives
    
    R_med{R_ctr,1} = get_nice_population_names(str2title(demographic_models_struct.model_str{i_d}));
    for j=1:length(show_s_null)
        [median_val median_ind] = min(abs( demographic_models_struct.data{i_d}.p_vec(j,:) - demographic_models_struct.data{1}.p_vec(j,end)/2) );
        demographic_models_struct.data{i_d}.median_vec(j) = demographic_models_struct.data{1}.x_vec(median_ind);
        R_med{R_ctr, length(show_s_null)-j+2} = num2str(demographic_models_struct.data{i_d}.median_vec(j), 2);
    end
    R_ctr=R_ctr+1;
end
my_mkdir(fullfile(new_figs_dir, 'all_models'));
savecellfile(R_med, fullfile(new_figs_dir, 'all_models', 'median_allele_all_models.txt'));

if(size(R_med,1) >= 13)
    R_med_inds = [1 7 9 6 8 11 13];
else
    R_med_inds = 1:4;
end
R_med_latex = latex(R_med(R_med_inds,1:end-1)', 2, precision); % Hard-coded populations indices
R_med_latex = mat2cell(R_med_latex, ones(size(R_med_latex,1),1));
R_med_latex{1} = strrep(R_med_latex{1},  '|c', '|r'); % align to right
R_med_latex{1} = strrep(R_med_latex{1},  '{|r', '{|l'); % align to left
savecellfile(R_med_latex, fullfile(new_figs_dir, 'all_models', 'median_allele_all_models_latex.txt'), [], 1);


populations_ordered = {'equil', 'expan1', 'expan2', '2phase', 'europ', 'finn1', 'finn2', 'ice', 'varsel1', 'varsel2', ...
    'finn1_bneck', 'finn1_nobneck', 'finn2_bneck', 'finn2_nobneck', 'ice_bneck', 'ice_nobneck'}; % Get ordering for plotting
[~, populations_perm, temp_perm] = intersect(demographic_models_struct.model_str, populations_ordered);
populations_perm = populations_perm(inv_perm(temp_perm));

if(length(demographic_models_struct.data) < demographic_models_struct.num_models)
    demographic_models_struct.data{demographic_models_struct.num_models} = [];
end
demographic_models_struct.file_names = vec2column(demographic_models_struct.file_names);
demographic_models_struct.model_str = vec2column(demographic_models_struct.model_str);
demographic_models_struct.data = vec2column(demographic_models_struct.data);
demographic_models_struct = order_struct_by_perm(demographic_models_struct, populations_perm, [], 2);  % new ordering
demographic_models_struct.num_models = length(populations_perm);

%save(equilibrium_parameters_output_file, '-append', 'demographic_models_struct'); % TEMP!

if(~figs_for_paper_flag) % plot many figures
    plot_frac_null_and_captured_internal(f_rare_vec, frac_null_by_freq_cumulative, xi_prob_leq_f, ...
        show_s_null_ind, rare_cumulative_per_gene, ...
        legend_vec, title_str, alpha_vec, figs_dir);
end % if ~figs_for_paper

close all;

n_samples = 5000; lambda=17.1; alpha=0.25; % NEW! Plot rho in finite sample ! (use LDLR values)
plot_rho_finite_size(demographic_models_struct.data{1}.x_vec, n_samples, lambda, alpha, ...
    demographic_models_struct, 5, new_figs_dir) % these give the population (Europe)


fig_files_to_ppt = plot_mean_num_alleles_and_het_vs_s_internal(demographic_models_struct, ...
    two_class_stat_struct, s_null_vec, N, mu_per_gene, show_s_null, ...
    color_vec, new_figs_dir);    % Plot mean number of alleles present in a single chromosome

for pop_groups = { {'equil'}, {'europ'}, ... % NEW: one main figure !!!
        {'equil', 'expan1', '2phase', 'europ', 'ice', 'finn1'}, ... % NEW: one figure for the paper for all 6 populations
        {'equil', 'expan1', 'expan2', '2phase'}, {'europ', 'ice', 'finn1', 'finn2'}, {'varsel1', 'varsel2'}, ...
        {'finn1_bneck', 'finn1_nobneck', 'finn2_bneck', 'finn2_nobneck'}, {'ice_bneck', 'ice_nobneck'} }; % new! add finland and iceland
    if(ismember('equil', pop_groups{1}))
        if(ismember('ice', pop_groups{1}))
            pop_str = 'all_models'; % all models for paper
        else
            if(length(pop_groups{1}) == 1)
                pop_str = 'equilibrium';
            else
                pop_str = 'expansion';
            end
        end
    else
        if(ismember('ice', pop_groups{1}))
            pop_str = 'bottleneck';
        else
            if(ismember('finn1_bneck', pop_groups{1}))
                pop_str = 'Finland';
            else
                if(ismember('ice_bneck', pop_groups{1}))
                    pop_str = 'Iceland';
                else
                    pop_str = 'varsel';
                end
            end
        end
        if(ismember('europ', pop_groups{1}) && (length(pop_groups{1}) == 1))
            pop_str = 'europe';
        end
    end
    tmp_fig_files_to_ppt = plot_frac_null_and_captured_by_population_internal(demographic_models_struct, pop_groups{1}, ...
        f_rare_vec, show_s_null, show_s_null_ind, w_x_null_mat, w_x_harmless, xi_prob_leq_f, my_symbol_vec, plot_bayes_factor, eric_color_vec, ...
        frac_null_by_freq_cumulative, pop_str, alpha_vec, figs_for_paper_flag, legend_vec, new_figs_dir);
    fig_files_to_ppt = union(fig_files_to_ppt, tmp_fig_files_to_ppt);
end % loop on pop. groups

close all;
plot_frac_null_and_harmless_all_populations_internal(demographic_models_struct, ...
    show_s_null_ind, xi_prob_leq_f, frac_null_by_freq_cumulative,  eric_color_vec, my_symbol_vec, legend_vec, title_str, ...,
    figs_for_paper_flag, figs_dir);
close all;

if(figs_for_paper_flag) % New: Plot cumulative variance explained
    plot_cumulative_var_explained_internal(show_s_null, f_rare_vec, N, figs_dir);
end % if plot figs for paper

%%%%%%%%%%%
% Make some simple plots of density
if(~figs_for_paper_flag)
    plot_simple_density_internal(f_rare_vec, s_null_vec, N, show_s_null_ind, legend_vec, new_figs_dir);
end % optional plots

close all;

if(~figs_for_paper_flag)
    % Plot isolines for variance explained with effect size and s on the two axis
    % %     GRR_vec = 1:0.01:8; % genetic relative risk
    f_vec = 0.0000001; % temp (should change f, but the resulting beta isn't very sensitive to allele freq.
    total_desired_var_explained = logspace(-4,-1, 7);
    plot_var_explained_by_s_and_grr_internal(demographic_models_struct, two_class_stat_struct, s_null_vec, N, mu, gene_length, ...
        total_desired_var_explained, f_vec, prevalence, color_vec, my_symbol_vec, new_figs_dir);
    
    
    plot_var_explained_vs_s_internal(two_class_stat_struct, s_null_vec, show_s_null, show_s_null_ind, N, alpha_vec, ...
        prob_null_given_x_vec, figs_dir, new_figs_dir); % New! Plot variance explained by one locus as a function of s. Also save table
    
    
    plot_power_proxy_internal(frac_null_by_freq_cumulative, show_s_null_ind, xi_prob_leq_f, f_rare_vec, alpha_vec, ...
        legend_vec, figs_dir);  % Plot power (what is proxy for power?)
    
    
end % if figs_for_paper

%%%%%%%%%%%%%%%%%%%%%% Figure 2: Plot things as function of s, the selection coefficient %%%%%%%%%%%%%%%%%%%%%%


close all;





% Save as powerpoint slides
%mxdom2ppt(fig_files_to_ppt,'ppt')
% publish(fig_files{1}, 'ppt')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tmp_fig_files_to_ppt = plot_frac_null_and_captured_internal(f_rare_vec, frac_null_by_freq_cumulative, xi_prob_leq_f, ...
    show_s_null_ind, rare_cumulative_per_gene, ...
    legend_vec, title_str, alpha_vec, figs_dir)
for log_x_flag = 0:1 % Plot frac null and frac captured on the same plot
    figure;
    tmp_title_str = [strrep(title_str, ...
        'Detect. power rare', 'Frac. null') ' ; frac. null (solid); frac. captured (dashed)'];
    tmp_title_str = strdiff(tmp_title_str, 'cum. freq.');
    %    title(tmp_title_str);
    switch log_x_flag
        case 0 % linear plot
            x_log_str = '';
            plot(f_rare_vec, frac_null_by_freq_cumulative{1}(show_s_null_ind,:), 'linewidth', 2); hold on;
            plot(f_rare_vec, xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene, '--', 'linewidth', 2);
        case 1 % logarithmic plot
            x_log_str = '_log';
            semilogx(f_rare_vec, frac_null_by_freq_cumulative{1}(show_s_null_ind,:), 'linewidth', 2); hold on;
            semilogx(f_rare_vec, xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene, '--', 'linewidth', 2);
    end
    
    legend(legend_vec, 4-log_x_flag);   legend('boxoff');
    %    legend('frac. null', 'frac mutations captured');
    xlabel('Derived allele frequency f'); ylabel('Frac. null (solid) ; Frac. captured (dashed)');
    ylim([0 1.01*max(max(frac_null_by_freq_cumulative{1}(show_s_null_ind,:)))]);
    my_saveas(gcf, fullfile(figs_dir, 'power', ...
        ['rare_variants_null_fraction_as_function_of_freq_cutoff_alpha_' num2str(alpha_vec,2) x_log_str]), {'epsc', 'pdf'});
end % loop on log flag
tmp_fig_files_to_ppt = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot fraction of nulls captured as function of threshold
function fig_files_to_ppt = plot_frac_null_and_captured_by_population_internal(demographic_models_struct, pop_group, ...
    f_rare_vec, show_s_null, show_s_null_ind, w_x_null_mat, w_x_harmless, xi_prob_leq_f, my_symbol_vec, plot_bayes_factor, eric_color_vec, ...
    frac_null_by_freq_cumulative, pop_str, alpha_vec, figs_for_paper_flag, legend_vec, new_figs_dir)

AssignGeneralConstants;

orange = [1 0.6 0.1];
population_color_vec = {'k', orange, 'r', 'g', 'b', 'c'};
eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
my_symbol_vec = {'--', '-'}; % flip ordering (set integer powers as solid lines)
selection_color_vec = {'k', 'b', 'g', orange, 'r'}; % replace yellow with orange

num_models_to_plot = length(pop_group);
num_plots_x = ceil(num_models_to_plot/2);
num_plots_y = min(num_models_to_plot, 2);


fig_files_to_ppt = [];
for plot_type = [0 2] % 0:0 % 1 % 0 - frac. which are null. 1 - frac. which are captured
    model_ctr=1;  figure;
    for i_d = 1:demographic_models_struct.num_models    % New!!! plot also for other populations !!!!
        neutral_ind = find(demographic_models_struct.data{i_d}.s_vec == 0)
        
        if(ismember(demographic_models_struct.model_str{i_d}, pop_group))
            num_s = size(demographic_models_struct.data{i_d}.p_vec, 1);
            
            [~, plot_s_inds] = intersect(demographic_models_struct.data{i_d}.s_vec, ...
                [10.^(-1:-0.5:-5) 0]); % plot_s_inds = plot_s_inds(end:-1:1);
            switch demographic_models_struct.model_str{i_d}
                case 'equilibrium'
                    plot_f_vec = f_rare_vec;
                otherwise
                    demographic_models_struct.data{i_d} = load(demographic_models_struct.file_names{i_d}); % this gives both s-vec and cumulatives
                    plot_f_vec = demographic_models_struct.data{i_d}.x_vec;
            end
            
            if(~strcmp(pop_str, 'all_models'))
                plot_model_ctr = model_ctr;
            else  % change horozontal to vertical order. Put xlabel only at bottom
                switch model_ctr
                    case 1
                        plot_model_ctr=1;
                    case 2
                        plot_model_ctr=3;
                    case 3
                        plot_model_ctr=5;
                    case 4
                        plot_model_ctr=2;
                    case 5
                        plot_model_ctr=4;
                    case 6
                        plot_model_ctr=6;
                end
            end
            
            switch plot_type
                case 0  % fraction which are null, rho
                    if(plot_bayes_factor)
                        y_vec = w_x_null_mat{1}(show_s_null_ind,:) ./ ...
                            repmat(w_x_harmless{1}, length(show_s_null_ind), 1); %   y_vec ./ alpha_vec;
                        y_str = 'relative depletion for nulls';
                    else
                        switch demographic_models_struct.model_str{i_d}
                            case 'equilibrium'
                                y_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:);
                            otherwise
                                if(num_s == 10)
                                    y_vec = alpha_vec .* demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) ./ ... % take null
                                        ( alpha_vec .* demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) + ... % divide by sum of nulls and neutrals
                                        (1-alpha_vec) .* repmat(demographic_models_struct.data{i_d}.p_vec(neutral_ind,:), num_s, 1) ); % This should go DOWN!!!
                                else % NEW! multiple
                                    y_vec = alpha_vec .* demographic_models_struct.data{i_d}.p_vec(1:end,:) .* repmat(demographic_models_struct.data{i_d}.num_alleles_per_chrom(1:end)', 1, length(plot_f_vec)) ./ ... % take null
                                        ( alpha_vec .* demographic_models_struct.data{i_d}.p_vec(1:end,:) .* repmat(demographic_models_struct.data{i_d}.num_alleles_per_chrom(1:end,:)', 1, length(plot_f_vec)) + ... % divide by sum of nulls and neutrals
                                        (1-alpha_vec) .* repmat(demographic_models_struct.data{i_d}.p_vec(neutral_ind,:).* demographic_models_struct.data{i_d}.num_alleles_per_chrom(neutral_ind) , num_s, 1) ); % This should go DOWN!!!
                                end
                                %                                 y_vec2 = frac_null_conditional_on_freq_less_f(  ...
                                %                                     0.001, [], alpha_vec, demographic_models_struct.data{i_d}.x_vec, ...
                                %                                     demographic_models_struct, demographic_models_struct.model_str{i_d}) % use new function
                                %                                 y_vec_finite_sample = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
                                %                                     0.001, [], alpha_vec, demographic_models_struct.data{i_d}.x_vec, ...
                                %                                     demographic_models_struct, demographic_models_struct.model_str{i_d}, 100) % use new function - finite sample
                        end
                        y_str = 'Proportion of null alleles, \rho_s(f)';
                    end
                    cur_file_str = 'rho_f'; % 'frac_null_given_freq_less_than_cutoff_alpha_'
                    cur_title_str = ['\rho(f) for \alpha=' num2str(alpha_vec)];
                    %                prob_null_given_x_vec = y_vec(:,end);
                case 1  %-cumulative fraction captured (union of nulls and neutrals)
                    switch demographic_models_struct.model_str{i_d}
                        case 'equilibrium'
                            y_vec = xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene;
                        otherwise
                            y_vec = bsxfun(@plus, demographic_models_struct.data{i_d}.p_vec(end:-1:1,:), ...
                                demographic_models_struct.data{i_d}.p_vec(end,:));
                    end % switch
                    cur_file_str = 'rare_variants_freq_less_than_cutoff_alpha_'
                    cur_title_str = ['Frac. of all alleles with freq. below f^*. \alpha=' num2str(alpha_vec)];
                    y_str = 'Proportion of alleles below f, \Psi_s(f)';
                    
                case 2 % cumulative fraction captured (only null)
                    y_str = 'Proportion of alleles below f, \Psi_s(f)';
                    switch demographic_models_struct.model_str{i_d}
                        case 'equilibrium'
                            plot_f_vec = f_rare_vec;
                            y_vec = cur_frac_vec .* ...
                                xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene;
                        otherwise
                            if((i_d>length(demographic_models_struct.data)) || (~isfield(demographic_models_struct.data{i_d}, 'x_vec')))
                                continue;
                            end
                            plot_f_vec = demographic_models_struct.data{i_d}.x_vec;
                            if(num_s == 10)
                                y_vec = demographic_models_struct.data{i_d}.p_vec; %   (end:-1:1,:); % FLIP !!!! who knows whats going on!!!
                            else
                                y_vec = demographic_models_struct.data{i_d}.p_vec; % FLIP !!!! who knows whats going on!!!
                            end
                    end
                    cur_file_str = 'cumul_f';
                    normalize_str = 'normalized_';
                    y_vec = y_vec ./ repmat(y_vec(:,end), 1, size(y_vec, 2)); % normalized to get to one
            end % switch plot type
            %%%%%%            cur_title_str = [cur_title_str ' . ' demographic_models_struct.model_str{i_d}];
            for log_x_flag = 1 % 0:1
                for log_y_flag = 0 % 0:1
                    if(num_models_to_plot > 1)
                        subplot(num_plots_x,2, plot_model_ctr); %  figure;     % New plot separately and on log scale
                    end
                    %%            if( (figs_for_paper_flag && log_x_flag) ) %  && (~log_y_flag) )
                    
                    %%            else
                    %%                subplot(2,2, log_x_flag + 2*log_y_flag+1);
                    %%            end
                    switch log_x_flag
                        case 0
                            x_log_str = '';
                            switch log_y_flag
                                case 0
                                    y_log_str = '';
                                    %    plot(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2); hold on;
                                    plot(f_rare_vec, y_vec, 'linewidth', 2);
                                case 1
                                    y_log_str = ' (log)';
                                    %    semilogy(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2);  hold on;
                                    semilogy(f_rare_vec,y_vec, 'linewidth', 2);
                            end
                        case 1
                            x_log_str = ' (log)';
                            switch log_y_flag
                                case 0
                                    y_log_str = '';
                                    %    semilogx(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2);  hold on;
                                    for j=1:length(plot_s_inds) %  size(y_vec, 1) % loop on all s values?
                                        semilogx(plot_f_vec, y_vec(plot_s_inds(j),:), ...
                                            'color', selection_color_vec{ceil(j/2)}, ...
                                            'linestyle', my_symbol_vec{mod_max(j,2)}, 'linewidth', 2); hold on;
                                        %                                            [eric_color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}],  'linewidth', 2); hold on;
                                    end
                                case 1
                                    y_log_str = ' (log)';
                                    %    loglog(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2);  hold on;
                                    for j=1:length(plot_s_inds) % size(y_vec, 1)
                                        loglog(plot_f_vec, y_vec(plot_s_inds(j),:), [eric_color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}], 'linewidth', 2); hold on;
                                    end
                            end
                    end % switch log x
                    if(~strcmp(pop_str, 'all_models') || ismember(plot_model_ctr, [5 6])) % put xlabel only at bottom
                        xlabel('Derived allele frequency f'); % x_log_str
                    end
                    if(~strcmp(pop_str, 'all_models') || ismember(plot_model_ctr, [3 4])) % put xlabel only at bottom
                        ylabel([y_str y_log_str]);
                    end
                    if(log_x_flag)
                        xlim([10^(-4) 1]);
                    end
                    if(log_y_flag) % set the same axis for x and y
                        ylim([10^(-4) 1]);
                    end
                    if(plot_type == 0)
                        ylim([0 alpha_vec .* 1.05]);
                    else
                        ylim([0 1.05]);
                    end
                    title(get_nice_population_names(str2title(demographic_models_struct.model_str{i_d}))); % no title !!
                    
                    %    ylim([0 1.01*max(max(frac_null_by_freq_cumulative(show_s_null_ind,:)))]);
                    if(figs_for_paper_flag)
                        switch plot_type
                            case 0
                                output_fig_file = 'figure1b';
                                title('Fig 1b', 'fontsize', 14, 'fontweight', 'b');
                                if(plot_bayes_factor)
                                    ylim([0 1.01]);
                                end
                            case 1 % plot combined null_neutral freq.
                                if(~log_y_flag)
                                    ylim([0 1]);
                                end
                                output_fig_file = ['figure1c' y_log_str];
                                title('Fig 1c', 'fontsize', 14, 'fontweight', 'b');
                        end
                    else
                        switch log_y_flag
                            case 0
                                if(plot_type == 0)
                                    output_fig_file = fullfile(demographic_models_struct.model_str{i_d}, ...
                                        ['fig2c_' cur_file_str num2str(alpha_vec,2)]);
                                else
                                    output_fig_file = fullfile(demographic_models_struct.model_str{i_d}, ...
                                        ['fig2b_' cur_file_str num2str(alpha_vec,2)]);
                                end
                                
                            case 1
                                output_fig_file = fullfile('power', [cur_file_str num2str(alpha_vec,2)]);
                        end % switch plot_type
                    end
                    %            if(~figs_for_paper_flag)
                    %            end
                    add_faint_grid(0.5);
                    
                    %%                    my_saveas(gcf, fullfile(new_figs_dir, output_fig_file), {'epsc', 'pdf'});
                end % loop on log-y
            end % loop on log-x
            model_ctr=model_ctr+1;
        end % if is-member
    end % loop on population models
    %        suptitle(cur_title_str); % make one title for all
    
    if(num_models_to_plot > 1)
        h_l = legend(legend_vec, 3, 'fontsize', 10); % 'location', 'east');
        pos_l = get(h_l, 'position'); set(h_l, 'position', [pos_l(1)+0.275 pos_l(2)-0.15 pos_l(3) pos_l(4)]);
        legend('boxoff'); % +plot_type*(~figs_for_paper_flag));
        orient landscape;
    end
    my_saveas(gcf, fullfile(new_figs_dir, ['all_models/' cur_file_str '_' pop_str]), {'epsc', 'pdf', 'jpg'});
end % loop on plot type


for log_y_flag=0:1   %New: plot quantiles
    model_ctr=1;  figure;
    %     quantiles_vec = [0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99];
    %     quantiles_legend_vec = num2str(quantiles_vec');
    
    if(log_y_flag)
        log_str = '_log';
    else
        log_str = '';
    end
    
    for i_d = 1:demographic_models_struct.num_models    % New!!! plot also for other populations !!!!
        
        if(ismember(demographic_models_struct.model_str{i_d}, pop_group))
            switch demographic_models_struct.model_str{i_d}
                case 'equilibrium'
                    plot_s_vec = s_null_vec;
                otherwise
                    demographic_models_struct.data{i_d} = load(demographic_models_struct.file_names{i_d}); % this gives both s-vec and cumulatives
                    f_null_vec = bsxfun(@rdivide, demographic_models_struct.data{i_d}.p_vec(end:-1:1,:), ...
                        demographic_models_struct.data{i_d}.p_vec(end:-1:1,end));  % make sure this is normalized. This should go UP!
                    plot_s_vec = show_s_null; % demographic_models_struct.data{i_d}.x_vec;
            end
            
            subplot(num_plots_x,2,model_ctr);
            quantiles_legend_vec = plot_quantiles(vec2row(demographic_models_struct.data{i_d}.x_vec), f_null_vec, plot_s_vec, log_y_flag);
            
            xlabel('Selection Coefficient s'); % x_log_str
            ylabel('Derived allele frequency f'); % x_log_str
            title(str2title(demographic_models_struct.model_str{i_d}));
            ylim([10^(-5) 1]);
            add_faint_grid(0.7); % add grid
            
            model_ctr = model_ctr+1;
            
        end % if is member
    end % loop on demographic models
    
    h_l = legend({quantiles_legend_vec}, 3, 'fontsize', 10); % 'location', 'east');
    pos_l = get(h_l, 'position'); set(h_l, 'position', [pos_l(1)+0.275 pos_l(2)-0.15 pos_l(3) pos_l(4)]); % change legend position
    legend('boxoff'); % +plot_type*(~figs_for_paper_flag));
    
    orient landscape;
    cur_file_str = 'quantiles_rare_variants_freq_less_than_cutoff_alpha_';
    my_saveas(gcf, fullfile(new_figs_dir, ['all_models/' cur_file_str log_str '_' pop_str]), {'epsc', 'pdf'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Plot SFS for different selection coefficients 
function plot_frac_null_and_harmless_all_populations_internal(demographic_models_struct, ...
    show_s_null_ind, xi_prob_leq_f, frac_null_by_freq_cumulative, eric_color_vec, my_symbol_vec, legend_vec, title_str, ...
    figs_for_paper_flag, figs_dir)
for i_d = 1:demographic_models_struct.num_models
    for log_x_flag = 1 % 0:1 % Plot null and harmless alleles separately
        for neutral_flag = 0 % 0:1 % 0: null; 1: neutral
            switch neutral_flag
                case 0
                    neutral_str = 'null';
                    cur_frac_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:);
                case 1
                    neutral_str = 'neutral';
                    cur_frac_vec = 1-frac_null_by_freq_cumulative{1}(show_s_null_ind,:);
            end
            switch demographic_models_struct.model_str{i_d}
                case 'equilibrium'
                    plot_f_vec = f_rare_vec;
                    cur_y_vec = cur_frac_vec .* ...
                        xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene;
                otherwise
                    if((i_d>length(demographic_models_struct.data)) || (~isfield(demographic_models_struct.data{i_d}, 'x_vec')))
                        continue;
                    end
                    plot_f_vec = demographic_models_struct.data{i_d}.x_vec;
                    cur_y_vec = demographic_models_struct.data{i_d}.p_vec(end:-1:1,:);
            end
            for normalize_flag = 1 % 0:1 % choose if to normalize (divide by
                switch normalize_flag
                    case 0 % do nothing
                        normalize_str = '';
                    case 1
                        normalize_str = 'normalized_';
                        cur_y_vec = cur_y_vec ./ repmat(cur_y_vec(:,end), 1, size(cur_y_vec, 2)); % normalized to get to one
                end
                figure;    % one plot
                tmp_title_str = [strrep(title_str, ...
                    'Detect. power rare', 'Frac. null captured') ' ;'];
                tmp_title_str = strdiff(tmp_title_str, 'cum. freq.');
                %    title(tmp_title_str);
                switch log_x_flag
                    case 0 % linear plot
                        x_log_str = '';
                        plot(plot_f_vec, cur_y_vec, 'linewidth', 2);
                        %             plot(f_rare_vec, frac_null_by_freq_cumulative{1}(show_s_null_ind,:), 'linewidth', 2); hold on;
                        %             plot(f_rare_vec, xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene, '--', 'linewidth', 2);
                    case 1 % logarithmic plot
                        x_log_str = '_log';
                        for j=1:size(cur_y_vec, 1)
                            semilogx(plot_f_vec, cur_y_vec(j,:), [eric_color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}], 'linewidth', 2); hold on;
                        end
                end
                %            legend_vec = [repmat('s^*= -', length(show_s_null), 1) num2str(show_s_null',3)];
                legend(legend_vec, 3-log_x_flag); legend('boxoff');
                %    legend('frac. null', 'frac mutations captured');
                xlabel('Derived allele frequency f'); ylabel('Proportion of alleles below f, \Psi_s(f)');
                ylim([0 1.01*max(max(cur_y_vec))]);
                if(log_x_flag)
                    xlim([10^(-4) 1]);
                end
                if(figs_for_paper_flag)
                    if(neutral_flag)
                        title('Fig 1b', 'fontsize', 14, 'fontweight', 'b');
                    else
                        ylabel('frac. alleles');
                        title('Fig 1a', 'fontsize', 14, 'fontweight', 'b');
                    end
                    output_fig_file = 'figure1a';
                else
                    %                 output_fig_file = fullfile('power', ['rare_variants_' normalize_str neutral_str '_fraction_captured_as_function_of_freq_cutoff_alpha_' ...
                    %                     num2str(alpha_vec,2) x_log_str]); % fig doesn't depend on alpha
                    output_fig_file = fullfile('new_eric', demographic_models_struct.model_str{i_d}, ...
                        ['fig2a_rare_variants_' normalize_str neutral_str '_cumulative_distribution_f_x' ...
                        x_log_str]);
                end
                add_faint_grid();
                title(tmp_title_str);
                my_saveas(gcf, fullfile(figs_dir, output_fig_file), {'epsc', 'pdf'});
            end % loop on normalization flag
        end % loop on neutral flag
    end % loop on log flag
end % loop on models




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_cumulative_var_explained_internal(show_s_null, f_rare_vec, N, figs_dir)
var_expl_vec = zeros(length(show_s_null), length(f_rare_vec));
for i=1:length(show_s_null)
    run_i = i
    for j=1:length(f_rare_vec)
        %    var_expl_vec(i,:) = XXX(f_rare_vec, g =
        %   var_expl_vec(i,:) =  exp(allele_freq_spectrum(f_rare_vec, -show_s_null(i), N, 0, 'log', 1));
        var_expl_vec(i,j)=  absorption_time_by_selection(show_s_null(i), 1, N, 0, f_rare_vec(j), -2);
        %        absorption_time_by_selection(s, theta, N, x_min, x_max, weight_flag)
    end
    var_expl_vec(i,:) = var_expl_vec(i,:) ./ var_expl_vec(i,end); % normalize
end
% for i=1:length(show_s_null)
%     var_expl_vec(i,:) =
figure;
%plot(f_rare_vec, var_expl_vec, 'linewidth', 2);
semilogx(f_rare_vec, var_expl_vec, 'linewidth', 2);
xlim([10^(-4) 1]);
xlabel('Derived allele frequency f'); ylabel('cum. var. expl.');
% legend(legend_vec, 3);
title('Fig 1d', 'fontsize', 14, 'fontweight', 'b');
my_saveas(gcf, fullfile(figs_dir, 'figure1d'), {'epsc', 'pdf'});























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function fig_files_to_ppt = plot_var_explained_by_s_and_grr_internal(demographic_models_struct, two_class_stat_struct, s_null_vec, N, mu, L, ...
    total_desired_var_explained, f_vec, prevalence, color_vec, my_symbol_vec, new_figs_dir)

fig_files_to_ppt = [];
for pop_groups = { {'equil', 'expan1', 'expan2', '2phase'}, {'europ', 'ice', 'finn1', 'finn2'}} % , {'equilibrium'} }
    for iso_plot_type = {'isocurves', 'many_lines'} % loop on how to plot var. expalined for different s and GRR
        if(ismember('equil', pop_groups{1}))
            pop_str = 'expansion';
        else
            if(ismember('equilibrium', pop_groups{1}))
                pop_str = 'equilibrium';
            else
                pop_str = 'bottleneck';
            end
        end
        model_ctr=1; figure;
        for i_d = 1:demographic_models_struct.num_models    % New!!! plot also for other populations !!!! Plot Isocurves
            switch demographic_models_struct.model_str{i_d}
                case 'equilibrium'
                    mean_het_vec_per_gene = 2.*N.*mu .* L .* 2 .* two_class_stat_struct.mean_het_x .* two_class_stat_struct.normalization_factor_x;
                    %                num_alleles_per_chr_vec = 2 .* N .* mu .* L .* 2 .* two_class_stat_struct.normalization_factor_x, 'b', 'linewidth', 2;
                    plot_s_null_vec = s_null_vec;
                    
                    % % %             case
                    % % %             for pop_groups = { {'equil', 'expan1', 'expan2', '2phase'}, {'europ', 'ice', 'finn1', 'finn2'} };
                    % % %     if(ismember('equil', pop_groups{1}))
                    % % %         pop_str = 'expansion';
                    % % %     else
                    % % %         pop_str = 'bottleneck';
                    % % %     end
                    
                    
                case {'equil', 'expan1', 'expan2', '2phase', 'europ', 'ice', 'finn1', 'finn2'}
                    plot_s_null_vec = demographic_models_struct.data{i_d}.s_vec; % only few val
                    demographic_models_struct.data{i_d}.density_p_vec = zeros(size(demographic_models_struct.data{i_d}.p_vec));
                    demographic_models_struct.data{i_d}.mean_het_x = zeros(length(plot_s_null_vec), 1);
                    for i_s = 1:length(plot_s_null_vec)
                        demographic_models_struct.data{i_d}.density_p_vec(i_s,:) = [vec2row(diff(demographic_models_struct.data{i_d}.p_vec(i_s,:))) 0];
                        demographic_models_struct.data{i_d}.mean_het_x(i_s) = ...
                            demographic_models_struct.data{i_d}.num_alleles_per_chrom(i_s) * ...
                            mean_hist( demographic_models_struct.data{i_d}.x_vec', demographic_models_struct.data{i_d}.density_p_vec(i_s,:));
                    end
                    mean_het_vec_per_gene = 2 .* (demographic_models_struct.data{i_d}.num_alleles_per_chrom - ...
                        demographic_models_struct.data{i_d}.mean_het_x); % demographic_models_struct.data{i_d}.num_alleles_per_chrom;
                otherwise
                    continue;
            end % switch demographic model
            if(~ismember(demographic_models_struct.model_str{i_d}, pop_groups{1}))
                continue;
            end
            
            GRR_isoline_vec = zeros(length(plot_s_null_vec), length(total_desired_var_explained));
            for j=1:length(total_desired_var_explained)
                desired_beta_vec = sqrt(total_desired_var_explained(j) ./ mean_het_vec_per_gene); % var. explained proportional to beta^2
                GRR_isoline_vec(:,j) = beta_to_genetic_relative_risk(desired_beta_vec, f_vec, prevalence);
            end
            
            figure; % subplot(2,2, model_ctr);
            switch iso_plot_type{1}
                case 'isocurves'
                    iso_plot_str = 'Isocurves for';
                    isocurve_legend_vec = cellstr( ...
                        [repmat('$V=10^{', length(total_desired_var_explained), 1) num2str(log10(total_desired_var_explained')) ...
                        repmat('}$', length(total_desired_var_explained), 1)])
                    for j=1:size(GRR_isoline_vec,2) % number of isolines
                        semilogx(plot_s_null_vec, GRR_isoline_vec(:,j), [color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}],  'linewidth', 2); hold on;
                    end
                    ylim([0.999 10]); ylabel('GRR');
                case 'many_lines'
                    iso_plot_str = '';
                    show_grr_vec = [1.3 1.5 2 3 5 10];
                    isocurve_legend_vec = cellstr([repmat('RR=', length(show_grr_vec), 1) num2str(show_grr_vec')]);
                    [~, show_beta_vec] = genetic_relative_risk_to_beta(f_vec, show_grr_vec, prevalence);
                    
                    if(strcmp(demographic_models_struct.model_str{i_d}, 'equil')) % plot equilibrium theoretical
                        mean_het_vec_per_gene = [vec2row(2.*N.*mu .* L .* 2 .* (2.*two_class_stat_struct.total_het_x)) ...
                            mean_het_vec_per_gene(1)]; %  .* two_class_stat_struct.normalization_factor_x;
                        plot_s_null_vec = [vec2row(s_null_vec) plot_s_null_vec(1)];
                        
                        no_nan_inds = find(~isnan(mean_het_vec_per_gene));
                        mean_het_vec_per_gene = mean_het_vec_per_gene(no_nan_inds);
                        plot_s_null_vec = plot_s_null_vec(no_nan_inds);
                        for j=1:length(show_beta_vec)
                            loglog(plot_s_null_vec, mean_het_vec_per_gene .* show_beta_vec(j).^2, color_vec(j), 'linewidth', 2); hold on; %  [color_vec(j) '--']
                        end
                    else % here plot Schaffner simulation
                        for j=1:length(show_beta_vec)
                            loglog(plot_s_null_vec, mean_het_vec_per_gene .* show_beta_vec(j).^2, color_vec(j), 'linewidth', 2); hold on;
                        end
                    end
                    
                    ylim([0 1]); xlim([10^(-5) 10^(-1)]); ylabel('Var. Explained');
            end % switch iso_plot_type
            xlabel('s');
            title(demographic_models_struct.model_str{i_d}, 'interpreter', 'latex');
            add_faint_grid(0.5);
            model_ctr=model_ctr+1;
            h_leg =  legend(isocurve_legend_vec, 3, 'interpreter', 'latex', 'fontsize', 10);
            set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
            my_saveas(gcf, fullfile(new_figs_dir, demographic_models_struct.model_str{i_d}, 'var_explained', ...
                [demographic_models_struct.model_str{i_d} '_var_explained_' iso_plot_type{1}]), ...
                {'epsc', 'fig', 'pdf'}); % put all figs in the same model
            
            switch pop_str
                case 'expansion'
                    fig_files_to_ppt{1} = fullfile(new_figs_dir, demographic_models_struct.model_str{i_d}, 'var_explained', ...
                        [demographic_models_struct.model_str{i_d} '_var_explained_' iso_plot_type{1}]);
            end
            
        end % loop on demographic models
        %         my_saveas(gcf, fullfile(new_figs_dir, demographic_models_struct.model_str{i_d}, ...
        %             'fig1b_var_explained_per_gene_as_function_of_s_and_GRR_isocurves'), {'epsc', 'pdf'}); % put all figs in the same model
        
        plot_two_by_two = 0;
        if(plot_two_by_two)
            h_l =  legend(isocurve_legend_vec, 3, 'interpreter', 'latex', 'fontsize', 10);
            pos_l = get(h_l, 'position'); set(h_l, 'position', [pos_l(1)+0.265 pos_l(2)-0.15 pos_l(3) pos_l(4)]);
            legend('boxoff'); % +plot_type*(~figs_for_paper_flag));
            suptitle([iso_plot_str ' variance explained per gene. \pi=' num2str(prevalence) '. ' pop_str]);
            orient landscape;
            my_saveas(gcf, fullfile(new_figs_dir, 'all_models/var_explained', ...
                [pop_str '_var_explained_' iso_plot_type{1}]), ...
                {'epsc', 'pdf'}); % put all figs in the same model
        end
    end % loop on iso plot type
end % loop on pop groups




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function plot_var_explained_vs_s_internal(two_class_stat_struct, s_null_vec, show_s_null, show_s_null_ind, N, alpha_vec, ...
    prob_null_given_x_vec, figs_dir, new_figs_dir) % New! Plot variance explained by one locus as a function of s
%tmp_num_loci_vec = -S_vec .* integral_phi_x' ./ (-S_vec - (1-exp(S_vec)))
tmp_var_explained_per_gene_vec = ... % This is the var explained per GENE
    allele_freq_cumulative(1-1/(2*N), -s_null_vec, N, 0, 'linear', 'var') - ...
    allele_freq_cumulative(1/(2*N), -s_null_vec, N, 0, 'linear', 'var'); % Start with variance explained per gene

tmp_var_explained_per_polymorphic_locus_vec = ... % This is the var explained per GENE
    ( allele_freq_cumulative(1-1/(2*N), -s_null_vec, N, 0, 'linear', 'var') - ...
    allele_freq_cumulative(1/(2*N), -s_null_vec, N, 0, 'linear', 'var') ) ./ ...
    ( allele_freq_cumulative(1-1/(2*N), -s_null_vec, N, 0, 'linear', 0) - ...
    allele_freq_cumulative(1/(2*N), -s_null_vec, N, 0, 'linear', 0) ); % Start with variance explained per gene
tmp_num_genes_vec = 1 ./ tmp_var_explained_per_gene_vec;
tmp_num_polymorphic_loci_vec = 1 ./ tmp_var_explained_per_polymorphic_locus_vec;


figure; loglog(s_null_vec, tmp_var_explained_per_gene_vec, 'linewidth', 2 ); hold on;
loglog(s_null_vec, tmp_var_explained_per_polymorphic_locus_vec, 'r', 'linewidth', 2 );
xlabel('-s'); ylabel('Var. Explained');
legend('gene', 'polymorphic allele');
my_saveas(gcf, fullfile(figs_dir, 'power', ...
    'var_explained_per_gene_as_function_of_selection'), ...
    {'epsc', 'pdf'});
figure; loglog(s_null_vec, tmp_num_genes_vec, 'linewidth', 2 ); hold on;
loglog(s_null_vec, tmp_num_polymorphic_loci_vec, 'r', 'linewidth', 2 );
xlabel('-s'); ylabel('Num. Genes');
legend({'gene', 'polymorphic allele'},4);
my_saveas(gcf, fullfile(figs_dir, 'power', ...
    'num_genes_needed_as_function_of_selection'), ...
    {'epsc', 'pdf'});

% figure; loglog(s_null_vec, tmp_var_explained_per_polymorphic_locus_vec, 'linewidth', 2 );
% xlabel('-s'); ylabel('Var. Explained');
% my_saveas(gcf, fullfile(figs_dir, 'power', ...
%     'rare_variants_var_explained_per_polymorphic_allele_as_function_of_selection'), ...
%     {'epsc', 'pdf'});
% figure; loglog(s_null_vec, tmp_num_polymorphic_loci_vec, 'linewidth', 2 );
% xlabel('-s'); ylabel('Num. Genes');
% my_saveas(gcf, fullfile(figs_dir, 'power', ...
%     'num_polymorphic_alleles_needed_as_function_of_selection'), ...
%     {'epsc', 'pdf'});



% Alternative: Numeric version
x_vec = (1:N-1) ./ (N);
s = -0.0001; y_vec = allele_freq_spectrum(x_vec, s, N, 0);
z_vec = normalize_hist(x_vec, y_vec);
VV = integral_hist(x_vec, x_vec .* (1-x_vec) .* z_vec), 1 / VV

two_class_stat_struct.median_x = max(two_class_stat_struct.median_x, 0); % can't have negative values
for plot_flag = 0:1 % what to plot
    for log_flag = 0:1
        figure;
        switch log_flag
            case 0
                log_str = '';
                switch plot_flag
                    case 0
                        plot(s_null_vec, two_class_stat_struct.mean_x, 'linewidth', 2); hold on;
                        plot(s_null_vec, two_class_stat_struct.median_x, 'r', 'linewidth', 2);
                    case 1
                        plot(s_null_vec, two_class_stat_struct.mean_het_x, 'g', 'linewidth', 2);
                end
            case 1
                log_str = '_log';
                switch plot_flag
                    case 0
                        semilogx(s_null_vec, two_class_stat_struct.mean_x, 'linewidth', 2); hold on;
                        semilogx(s_null_vec, two_class_stat_struct.median_x, 'r', 'linewidth', 2);
                    case 1
                        semilogx(s_null_vec, two_class_stat_struct.mean_het_x, 'g', 'linewidth', 2);
                end
        end
        %title('Mean allele freq. as function of s');
        switch plot_flag
            case 0
                het_str = 'allele_freq';
                legend('mean-allele-freq.', 'median-allele-freq.');
            case 1
                het_str = 'heterozygosity';
                legend('mean-heterozygosity');
        end
        xlabel('-s'); % ylabel('mean-allele-freq.');
        my_saveas(gcf, fullfile(new_figs_dir, 'power', ...
            ['mean_' het_str '_given_s_' num2str(alpha_vec,2) log_str] ), {'epsc', 'pdf'});
    end % loop on log-flag
end % loop on plot flag

R = [ show_s_null'  vec2column(two_class_stat_struct.mean_x(show_s_null_ind)) ...
    vec2column(two_class_stat_struct.median_x(show_s_null_ind)) vec2column(two_class_stat_struct.mean_het_x(show_s_null_ind)) prob_null_given_x_vec]; % Save also a table
R = [{'S', 'mean-allele-freq', 'median-allele-freq.', 'mean-heterozygosity', ['prob-null (\alpha=' num2str(alpha_vec) ')']}' num2str_cell(num2cell(R'),2)]';
savecellfile(R, fullfile(new_figs_dir, 'power', 'mean_allele_freq_tab.txt'));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the total number of derived alleles per chromosome, CAF
function fig_files_to_ppt = plot_mean_num_alleles_and_het_vs_s_internal(demographic_models_struct, ...
    two_class_stat_struct, s_null_vec, N, mu_per_gene, show_s_null, ...
    color_vec, new_figs_dir)    % Plot mean number of alleles present in a single chromosome

AssignGeneralConstants; AssignRVASConstants;
mu=mu_per_site; L=625; %% mu=1.6*10^(-8); L=1500; % take typical parameters for mutation rate and gene length. We want to get mu*L=10^(-5)

for log_y_flag=0:1 % loop on log-flag
    switch log_y_flag
        case 0
            log_str = '';
        case 1
            log_str = '_log';
    end
    for mean_num_fig = [1 3]  % 2 % plot mean # alleles and mean heterozygosity as function of s
        switch mean_num_fig
            case 1
                mean_save_str = 'num_NULL_alleles'; allele_type_str = 'NULL';
            case 2
                mean_save_str = 'heterozygosity'; allele_type_str = 'HET';
            case 3
                mean_save_str = 'num_LOF_alleles'; allele_type_str = 'LOF';
        end
        for plot_curve = {'total_allele_freq', 'median_allele_freq'}
            figure; % New! plot in this figure all populations together
            switch plot_curve{1}
                case 'total_allele_freq'
                    symbol_str = '-';
                case 'median_allele_freq'
                    symbol_str = '-';
            end
            
            models_legend_vec = []; model_ctr=1;
            for i_d = 1:demographic_models_struct.num_models    % New!!! plot also for other populations !!!!
                if(ismember(demographic_models_struct.model_str{i_d}, ...
                        {'equil', 'europ', 'expan1', '2phase', 'finn1', 'ice'})) % 'expan2', 'finn2',
                    switch demographic_models_struct.model_str{i_d}
                        case 'equilibrium'
                            num_alleles_per_chr_vec = 2 .* N .* mu .* L .* 2 .* ...
                                two_class_stat_struct.normalization_factor_x; % , 'b', 'linewidth', 2;
                            plot_s_null_vec = s_null_vec;
                        case {'varsel1', 'varsel2', 'varselection1', 'varselection2'}
                            continue; % we don't plot for variable selection
                        otherwise
                            plot_s_null_vec = max(10^(-6), demographic_models_struct.data{i_d}.s_vec); % add a point at 10^(-6) to represent zero
                            switch plot_curve{1}
                                case 'total_allele_freq'
                                    num_alleles_per_chr_vec = demographic_models_struct.data{i_d}.num_alleles_per_chrom;
                                case 'median_allele_freq'
                                    num_alleles_per_chr_vec = demographic_models_struct.data{i_d}.median_vec;  % What does this mean??
                                    
                                    
                            end
                    end % switch demographic model
                    %                    symbol_str = '-';
                    %                     switch demographic_models_struct.model_str{i_d}
                    %                         case {'equil', 'expan1', 'expan2', '2phase'}
                    %                             symbol_str = '--';
                    %                         case {'europ', 'ice', 'finn1', 'finn2'}
                    %                             symbol_str = '-';
                    %                     end
                    
                    if(mean_num_fig == 3) % LOF alleles. Change overall rate (but same shape)
                        num_alleles_per_chr_vec = num_alleles_per_chr_vec .* (mu_per_gene.LOF / mu_per_gene.NULL);
                        rate_per_gene_str = num2str(mu*L .* (mu_per_gene.LOF / mu_per_gene.NULL), 3);
                    else
                        rate_per_gene_str = num2str(mu*L, 3);
                    end
                    switch log_y_flag
                        case 0
                            semilogx(plot_s_null_vec, num_alleles_per_chr_vec,  [color_vec(model_ctr) symbol_str], 'linewidth', 2); hold on;
                        case 1
                            loglog(plot_s_null_vec, num_alleles_per_chr_vec,  [color_vec(model_ctr) symbol_str], 'linewidth', 2); hold on;
                    end
                    models_legend_vec{model_ctr} = get_nice_population_names(demographic_models_struct.model_str{i_d});
                    num_alleles_per_chr_mat(model_ctr,:) = num_alleles_per_chr_vec;
                    model_ctr=model_ctr+1; % update model ctr
                end % if member of certain models
            end % loop on demographic model
            
            
            xlabel('Selection coefficient s', 'fontsize', 15, 'fontweight', 'bold');
            
            switch plot_curve{1}
                case 'total_allele_freq'
                    ylabel('Combined allele frequency f_{s}', 'fontsize', 15, 'fontweight', 'bold'); % Expected # of derived alleles
                    %                    title(['Expected ' str2title(mean_save_str) '  derived alleles per chrom. $\mu_{G,' allele_type_str '}=' rate_per_gene_str '$. '], ...
                    %                        'interpreter', 'latex', 'fontsize', 15);
                    fig_files_to_ppt{1} = fullfile(new_figs_dir, 'all_models', ...
                        ['fig1a_mean_' mean_save_str '_per_chromosome_as_function_of_s' log_str]); % save to ppt
                case 'median_allele_freq'
                    ylabel('Median Allele Frequency');
                    title(['Median Allele Frequency of randomly sampled alleles '], 'interpreter', 'latex');
                    fig_files_to_ppt{1} = fullfile(new_figs_dir, 'all_models', ...
                        ['fig1b_median_allele_freq_as_function_of_s' log_str]); % save to ppt
            end % switch
            % New! also plot the median allele frequency
            
            set(gca, 'fontsize', 15);
            add_faint_grid(0.5);
            
            h_leg = legend(models_legend_vec); % legend('boxoff');
            set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
            my_saveas(gcf, fig_files_to_ppt{1}, {'epsc', 'fig', 'pdf'});
            
            % Prepare table in text
            num_populations = 6; num_selection = 10;
            R = cell(num_populations+1, num_selection + 1);
            
            for k=1:num_selection
                R{1,k+1} = num2str(s_null_vec(k), 2);
            end
            for j=1:num_populations
                R{j,1} = models_legend_vec{j};
                for k=1:num_selection
                    R{j,k+1} = num2str(num_alleles_per_chr_mat(j,k), 2);
                end
            end
            
            S_cell = [];
            S_cell{1} = 's';
            for j=1:length(show_s_null)
                S_cell{j+1} = ['$10^{' num2str(log10(show_s_null(j)),3) '}$'];
            end
            S_cell{2} = '$0$'; % manual - neutral alleles
            R = [S_cell([1 end:-1:2])' R']';
            savecellfile(R, fullfile(new_figs_dir, 'all_models', [plot_curve{1} '_tab.txt'])); % already have
            
            R_CAF_latex = latex(R(1:7,[1 end:-1:2])', 2, precision);
            R_CAF_latex = mat2cell(R_CAF_latex, ones(size(R_CAF_latex,1),1));
            R_CAF_latex{1} = strrep(R_CAF_latex{1},  '|c', '|r'); % align to right
            R_CAF_latex{1} = strrep(R_CAF_latex{1},  '{|r', '{|l'); % align to left
            savecellfile(R_CAF_latex, fullfile(new_figs_dir, 'all_models', [plot_curve{1} '_tab_latex.txt']), [], 1);
            %        'combined_allele_freq_all_models_latex.txt'), [], 1);
            
        end % loop on plot curve
        
    end % loop on plot type
end % loop on log flag



% R_CAF = cell(demographic_models_struct.num_models+2, length(show_s_null)+2); % New! prepare table
% R_CAF{1,1} = 's';
% R_CAF{3,1} = 'f_{null}(equil.)';
% for j=1:length(show_s_null)
%     R_CAF{1, j+1} = ['$10^{' num2str(log10(show_s_null(j)),3) '}$'];
%     R_CAF{3, j+1} = num2str(n_null_vec_equilibrium(j), 3);
% end
% R_CAF{1,2} = '0';
%
% R_CAF{5,1} = 'Combined Allele Freq.:'; R_ctr=6;
% for i_d = 1:demographic_models_struct.num_models
%     run_id = i_d
%     if(i_d > length(demographic_models_struct.file_names))
%         continue;
%     end
%     if(~ismember(demographic_models_struct.model_str{i_d}, {'equil', '2phase', 'expan1', 'expan2', 'europ', 'ice', 'finn1', 'finn2'}))
%         continue;
%     end
%     demographic_models_struct.data{i_d} = load(demographic_models_struct.file_names{i_d}); % this gives both s-vec and cumulatives
%
%     R_CAF{R_ctr,1} = get_nice_population_names(str2title(demographic_models_struct.model_str{i_d}));
%     for j=1:length(show_s_null)
%         [median_val median_ind] = min(abs( demographic_models_struct.data{i_d}.p_vec(j,:) - demographic_models_struct.data{1}.p_vec(j,end)/2) );
%         demographic_models_struct.data{i_d}.median_vec(j) = demographic_models_struct.data{1}.x_vec(median_ind);
%         R_med{R_ctr, length(show_s_null)-j+2} = num2str(demographic_models_struct.data{i_d}.median_vec(j), 2);
%     end
%     R_ctr=R_ctr+1;
% end
% my_mkdir(fullfile(new_figs_dir, 'all_models'));
% savecellfile(R_CAF, fullfile(new_figs_dir, 'all_models', 'combined_allele_freq_all_models.txt'));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_power_proxy_internal(frac_null_by_freq_cumulative, show_s_null_ind, xi_prob_leq_f, f_rare_vec, alpha_vec, ...
    legend_vec, figs_dir)  % Plot power (what is proxy for power?)
y_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:).^2 .* xi_prob_leq_f{1}(show_s_null_ind,:);
y_vec = y_vec ./ repmat(max(y_vec,[],2), 1, size(y_vec, 2));
figure; plot(f_rare_vec, y_vec, 'linewidth', 2);
legend(legend_vec,4); title('Power proxy = (frac. null.)^2 * frac-kept, as function of allele-freq. cutoff');
xlabel('Derived allele frequency f'); ylabel('Proxy for power');
my_saveas(gcf, fullfile(figs_dir, 'power', ...
    ['rare_variants_power_vs_freq_cutoff_alpha_'  num2str(alpha_vec,2)]), {'epsc', 'pdf'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Plot Function 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_simple_density_internal(f_rare_vec, s_null_vec, N, show_s_null_ind, legend_vec, new_figs_dir)
cumulative_psi_s = cell(3,1);
f_rare_vec(end) = 0.999999999999;
for plot_type=2 % loop on different
    ctr=1;
    for j=show_s_null_ind
        cumulative_psi_s{plot_type}(ctr,:) = phi_s_integral(f_rare_vec,-s_null_vec(j)*4*N, plot_type-1);
        cumulative_psi_s{plot_type}(ctr,:) = ...
            cumulative_psi_s{plot_type}(ctr,:) - cumulative_psi_s{plot_type}(ctr,1);
        cumulative_psi_s{plot_type}(ctr,:) = ...
            cumulative_psi_s{plot_type}(ctr,:) ./ cumulative_psi_s{plot_type}(ctr,end); % normalize to one
        ctr=ctr+1;
        switch plot_type
            case 1 % here do \int_{0}^x psi_s(x)
                plot_str = '\phi_s(x)';
                
            case 2 % here do \int_0^x x*psi_s(x)
                plot_str = 'x\phi_s(x)';
                
            case 3 % here do \int_0^x x^2*psi_s(x)
                plot_str = 'x^2\phi_s(x)';
                
        end % switch
    end % loop on selection coefficient s
    file_str = strrep(plot_str, ' ', '_');
    file_str = strrep(file_str, '\', '_');
    file_str = strrep(file_str, '^', '_');
    file_str = strrep(file_str, '(', '_');
    file_str = strrep(file_str, ')', '_');
    
    figure; % Plot ..
    for log_x = 0:1
        for log_y = 0:1
            log_y_str = '';
            subplot(2,2,log_x+log_y*2+1);
            switch log_x
                case 0
                    log_x_str = '';
                    switch log_y
                        case 0
                            plot(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                            
                        case 1
                            log_y_str = '(log)';
                            semilogy(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                    end
                case 1
                    log_x_str = '(log)';
                    switch log_y
                        case 0
                            log_x_str = '(log)';
                            semilogx(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                            
                        case 1
                            log_y_str = '(log)';
                            loglog(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                            legend(legend_vec, 4);
                    end
            end
            title(['Cumulative fraction of distribution $\int_{0}^{f^*} ' plot_str ' dx$'], 'interpreter', 'latex');
            xlabel(['Derived allele Freq. f^* ' log_x_str]); ylabel(['Cumulative density ' log_y_str]);
        end % log x
    end % log y
    my_saveas(gcf, fullfile(new_figs_dir, 'power', ...
        ['cumulative_allele_freq_' file_str]), {'epsc', 'pdf'});
end % loop on plot type



