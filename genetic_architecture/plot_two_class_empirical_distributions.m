% Plot distribution of # allele, not just means
%
function plot_two_class_empirical_distributions(two_class_data_dir, allele_freq_file, figs_dir)

AssignGeneralConstants;
orange = [1 0.6 0.1];
population_color_vec = {'k', orange, 'r', 'g', 'b', 'c'};
eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
my_symbol_vec = {'--', '-'}; % flip ordering (set integer powers as solid lines)
selection_color_vec = {'k', 'b', 'g', orange, 'r'}; % replace yellow with orange


x_vec = ReadDataFile(allele_freq_file); x_vec = x_vec.derived_freq; % Note: allele freq. is in logarithmic coordinates.
% s_vec =  [10.^(-1:-0.5:-5) 0]; % , s=10^-1.5, s=10^-2, ... s=10^-5, s=0); % hard-coded selection based on Steve's mail. (discrapancy: 11 or 10 s?)

TwoClassFiles = GetFileNames(fullfile(two_class_data_dir, '*.mat'));
TwoClassFiles = setdiff(TwoClassFiles, GetFileNames(fullfile(two_class_data_dir, '*age*.mat')));
num_bins = 500; smooth_width = 5;

s_vec = zeros(1, length(TwoClassFiles)); % initilize all vector variables (set space)
mean_vec = s_vec; num_iters = s_vec;
n_vec = cell(1, length(TwoClassFiles)); p_vec = n_vec; n_hist_vec = n_vec; pop_vec = n_vec; median_vec = n_vec; f_null_vec = n_vec;
for i=1:length(TwoClassFiles) % go over each file
    load_file = i
    load(fullfile(two_class_data_dir, TwoClassFiles{i}));
    s_vec(i) = 10^-str2num( [str2word('.', TwoClassFiles{i}, 2) '.' str2word('.', TwoClassFiles{i}, 3)] );
    pop_vec{i} = str2word('.', TwoClassFiles{i}, 1);
    n_vec{i} = sum(s_distribution'); % just get total cumulative allele frequency
    [p_vec{i}, n_hist_vec{i}] = hist_density(n_vec{i}, num_bins, [], [], [], smooth_width);
    mean_vec(i) = mean(sum(s_distribution')); std_vec(i) = std(sum(s_distribution'));
    num_iters(i) = length(n_vec{i});
    
    s_cumulative = cumsum(s_distribution, 2);
    f_null_vec{i} = mean(s_cumulative); f_null_vec{i} = f_null_vec{i} ./ f_null_vec{i}(end);
    
    median_vec{i} = zeros(num_iters(i), 1);
    for j=1:num_iters(i)
        median_ind = find(s_cumulative(j,:) >= 0.5*s_cumulative(j,end), 1);
        median_vec{i}(j) = x_vec(median_ind); % NEW! Compute also median allele frequency !!!
    end
    
    median_vec_non_zero(i) = x_vec( find(f_null_vec{i}>=0.5, 1) ); % NEW! take median of just non-zero frequencies
end
s_vec(s_vec == 1) = 0; % format exception: 0 is denoted 0.0.

my_mkdir(fullfile(two_class_data_dir, 'all'));
save(fullfile(two_class_data_dir, 'all', 'all_models_statistics.mat'), ...
    's_vec', 'pop_vec', 'n_vec', 'p_vec', 'mean_vec', 'std_vec', 'num_iters', 'f_null_vec', 'median_vec');

six_populations_vec = {'equil', 'expan1', '2phase', 'europ', 'finn3', 'ice3'};
for population_str = six_populations_vec % {'equil', 'expan1', '2phase', 'europ', 'finn1', 'ice'} % , 'finn3', 'ice3'} % NEW! Save one file per population !!!!
    model_inds = strmatch(population_str{1}, TwoClassFiles);
    demographic_models_struct = [];
    demographic_models_struct.x_vec = x_vec;
    demographic_models_struct.s_vec = s_vec(model_inds);
    demographic_models_struct.p_vec = zeros(length(model_inds), length(x_vec));
    i_ctr=1;
    if(isempty(model_inds))
        continue;
    end
    for i=vec2row(model_inds)
        demographic_models_struct.p_vec(i_ctr,:) = f_null_vec{i};
        demographic_models_struct.num_alleles_per_chrom(i_ctr) = mean(n_vec{i}); % cumulative allele freq.
        median_ind = find(f_null_vec{i} >= 0.5*f_null_vec{i}(end), 1);
        demographic_models_struct.median_vec(i_ctr) = x_vec(median_ind); % NEW! Compute also median allele frequency !!!
        i_ctr=i_ctr+1;
    end
    save(fullfile(two_class_data_dir, 'all', ['cumul_' population_str{1} '_summary.mat']), ...
        '-struct', 'demographic_models_struct');
end

w_smooth = 4; % smoothing parameter
show_s_vals = [0 logspace(-5, -1, 9)]; % [0 logspace(-3, -1, 5)]; % New: plot only some values of s % logspace(-4, -1, 7)
for figure_type = 1:4 % New: plot median values just from the data
    R_freq = []; % structure to save to a cell file
    for j=1:length(show_s_vals)
        R_freq{j+1, 1} = ['$10^{' num2str(log10(show_s_vals(j)),3) '}$'];
    end
    R_freq{1,2} = '0';
    R_freq{1,1} = 's';
    figure;
    pop_ctr=1;
    for population_str = six_populations_vec %{'equil', 'expan1', '2phase', 'europ', 'finn1', 'ice'} % 'finn3', 'ice3'} % , 'finn3', 'ice3'}
        pop_legend_vec{pop_ctr} = get_nice_population_names(population_str{1});
        R_freq{1,pop_ctr+1} = pop_legend_vec{pop_ctr};
        model_inds = strmatch(population_str{1}, TwoClassFiles);
        cur_s_vec = max(10^(-6), s_vec(model_inds));
        [cur_s_vec sort_perm] = sort(cur_s_vec);
        [intersect_s s_I s_J] = intersect(show_s_vals, s_vec(model_inds));
        
        switch figure_type
            case 1
                save_freq = mean_vec(model_inds);
                loglog(cur_s_vec, smooth(mean_vec(model_inds(sort_perm)), w_smooth), ...
                    'color', population_color_vec{pop_ctr}, 'linewidth', 2); hold on;
            case 2
                loglog(cur_s_vec, smooth(median_cell(median_vec(model_inds(sort_perm))), w_smooth), ...
                    'color', population_color_vec{pop_ctr}, 'linewidth', 2); hold on;
            case 3
                loglog(cur_s_vec, smooth(mean_cell(median_vec(model_inds(sort_perm))), w_smooth), ...
                    'color', population_color_vec{pop_ctr}, 'linewidth', 2); hold on;
                save_freq = median_vec_non_zero(model_inds);
            case 4 % NEW! Take median IAF conditional on being > 0 (i.e. only for polymorphic!)
                save_freq = median_vec_non_zero(model_inds);
                loglog(cur_s_vec, smooth(median_vec_non_zero(model_inds(sort_perm)), w_smooth), ...
                    'color', population_color_vec{pop_ctr}, 'linewidth', 2); hold on;                
        end % switch figure type
        for j=1:length(show_s_vals)
            R_freq{j+1, pop_ctr+1} = num2str(save_freq(s_J(j)), 3);
        end
        
        pop_ctr=pop_ctr+1;
    end
    
    xlabel('Selection coefficient s', 'fontsize', 14);
    
    switch figure_type % ALWAYS add legend !!!
        case 1 % CAF
            ylabel('Combined allele frequency f_s', 'fontsize', 14);
            save_file = 'CAF_different_populations';
        case 2 % median of median allele frequency
            ylabel('Median allele frequency', 'fontsize', 14);
            save_file = 'Median_IAF_different_populations';
        case 3 % mean of median allele frequency
            ylabel('Mean of median allele frequency f_s', 'fontsize', 14);
            save_file = 'MeanMedian_IAF_different_populations';
        case 4
            ylabel('Median allele frequency f_s', 'fontsize', 14);
            save_file = 'Median_nonzero_IAF_different_populations';
    end % switch figure type
    h_leg = legend(pop_legend_vec, 3); % 'eastoutside'); % legend('boxoff'); % just legend
    set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
    ylim([10^(-6) 1]);
    add_faint_grid(0.5);
    my_saveas(gcf, fullfile(figs_dir, save_file), {'epsc', 'pdf', 'jpg'}); % NEW: add .jpg for Robert
    
    R_freq_latex = latex(R_freq(:,:), 2, precision); % Hard-coded populations indices
    R_freq_latex = mat2cell(R_freq_latex, ones(size(R_freq_latex,1),1));
    R_freq_latex{1} = strrep(R_freq_latex{1},  '|c', '|r'); % align to right
    R_freq_latex{1} = strrep(R_freq_latex{1},  '{|r', '{|l'); % align to left
        
    savecellfile(R_freq', [fullfile(figs_dir, save_file) '.txt'], [], 1);
    savecellfile(R_freq_latex, [fullfile(figs_dir, save_file) '_latex.txt'], [], 1);
end % loop on figure type


% Get rid of variable selection models
bad_inds = strmatch('varsel', pop_vec); bad_inds = union(bad_inds, strmatch('just', pop_vec));
good_inds = setdiff(1:length(TwoClassFiles), bad_inds);
s_vec = s_vec(good_inds); pop_vec = pop_vec(good_inds);
p_vec = p_vec(good_inds); n_hist_vec = n_hist_vec(good_inds); mean_vec = mean_vec(good_inds);
num_files = length(good_inds);


s_vec_str = cellstr([repmat('10^{', num_files, 1) num2str(log10(vec2column(s_vec)), 2) repmat('}', num_files, 1)]);
for i=find(s_vec == 0)
    s_vec_str{i} = '0';
end
s_vec_str = strrep_cell(s_vec_str, ' ', '');

for i=1:num_files
    legend_vec{i} = [get_nice_population_names(pop_vec{i})]; legend_vec{i} = [legend_vec{i} repmat(' ', 1, 6-length(legend_vec{i})) ' s=' s_vec_str{i}];
    legend_vec{i} = [legend_vec{i} repmat(' ', 1, 17-length(legend_vec{i}))];
    legend_vec{i} = [legend_vec{i} ' \mu=' num2str(mean_vec(i), 3)];
    legend_vec{i} = [legend_vec{i} repmat(' ', 1, 30-length(legend_vec{i}))];
    legend_vec{i} = [legend_vec{i} ' \sigma=' num2str(std_vec(i), 3)];
    s_ind = strfind(legend_vec{i}, 's=')
    mu_ind = strfind(legend_vec{i}, '\mu');
    %     legend_vec{i} =
    
    legend_vec_no_s{i} = [legend_vec{i}(1:s_ind-1) ' ' legend_vec{i}(mu_ind:end)];
    legend_vec_no_s{i} = strrep(legend_vec_no_s{i}, 'n1', 'n'); % don't need 1
    legend_vec_no_pop{i} = legend_vec{i}(s_ind:end);
    legend_vec_only_s{i} = str2word(' ', legend_vec_no_pop{i}, 1)
end
unique_pop_vec = unique(pop_vec);
%unique_pop_vec = unique_pop_vec([2 4 1 3 5 6 ]); % For coefficient of variation: Equilibrium,  Expan1, Expan2, Europe, Finland, Iceland % temp!!!! [2 3 4 6 1 8 ]
unique_pop_vec = unique_pop_vec([2 3 4 5 1 6 ]); % For other figures: Equilibrium, Europe, Expan1, Finland, Expan2, Iceland % temp!!!! [2 3 4 6 1 8 ]
num_models = length(unique_pop_vec);
unique_s_vec = unique(s_vec);

mu_factor = 1; % 0.17; % TEMP! MOVE from 10^(-5) to 1.7*10^(-6) for diruptive alleles

SI_flag = 1; % plot for SI
figure;
main_selection_color_vec = selection_color_vec;
if(~SI_flag)
    main_selection_color_vec{3} = 'r';
    main_selection_color_vec{5} = [1 0.8 0.05]; % light orange/yellow
end
for i_d = 1:num_models % Loop on models and plot separately for each population
    model_inds = strmatch(unique_pop_vec{i_d}, pop_vec);
    model_str = get_nice_population_names(unique_pop_vec{i_d});
    
    
    R_model = []; s_null_vec = [0 logspace(-5, -1, 9)];
    R_model{1,1} = 's'; R_model{2,1} = '0';
    for j=1:9
        R_model{j+2,1} = ['10^{-' num2str((11-j)/2) '}'];
    end
    s_null_vec_inds = intersect(model_inds, find(ismember(s_vec, s_null_vec)));
    
    quantiles_vec = [0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99];
    R_model{1,2} = 'Mean'; % R_model{1,3} = 'Median';
    R_model{1,3} = 'St.d.';
    for j=1:length(quantiles_vec) % get quantiles
        R_model{1,j+3} = ['$' num2str(100*quantiles_vec(j)) '\%$'];
    end
    
    for j=1:10
        R_model{j+1,2} = num2str(mu_factor .* mean_vec(s_null_vec_inds(j)), 2); % get mean
        %        R_model{j+1,3} = num2str(median_vec(s_null_vec_inds(j)), 2); % get median
        R_model{j+1,3} = num2str(sqrt(mu_factor) .* std_vec(s_null_vec_inds(j)), 2); % get std
        for k=1:length(quantiles_vec)
            R_model{j+1,k+3} = num2str( mu_factor .* quantile(n_vec{s_null_vec_inds(j)}, quantiles_vec(k)), 2);
        end
    end
    R_model = [R_model(:,1) R_model([1 2 11:-1:3] ,2:end)];
    my_mkdir(fullfile(figs_dir, 'different_models/variation'));
    savecellfile(R_model, fullfile(figs_dir, 'different_models/variation', ...
        [unique_pop_vec{i_d} '_quantiles.txt'] ));
    
    R_model_latex = latex(R_model', 2, 2);
    R_model_latex = mat2cell(R_model_latex, ones(size(R_model_latex,1),1));
    R_model_latex{1} = strrep(R_model_latex{1},  '|c', '|r'); % align to right
    R_model_latex{1} = strrep(R_model_latex{1},  '{|r', '{|l'); % align to left
    %    R_model_latex{2} = strrep(R_model_latex{2}, '%', '\%');
    R_model_latex{2} = strrep(strrep(strrep(strrep(R_model_latex{2}, '10^', '$10^'), '(', '{'), '}', '}$'), '10^{-6}', '0');
    savecellfile(R_model_latex, fullfile(figs_dir, 'different_models/variation', ...
        [unique_pop_vec{i_d} '_quantiles_tab_latex.txt'] ), [], 1);
    
    
    % Selection coefficient
    
    for plot_type = {'inverse_cumulative_normalized'} % inverse_cumulative_normalized'} % 'CAF_cumulative_log', 'inverse_cumulative_normalized', 's_confidence_interval', 'CAF_cumulative_log'} % , 'alpha_confidence_interval'} % , 'density', 'median_allele_frequency', 's_estimator_variance'}
        switch plot_type{1}
            case 'CAF_cumulative_log'
                subplot(3,2,i_d);
            case 'coefficient_of_variation' % do nothing (one curve per population) 
                
            case 'inverse_cumulative_normalized' % TEMP: add !!!
                if(SI_flag)
                    subplot(3,2,i_d);
                else
                    figure;
                end
            otherwise
                figure;
        end
        if(SI_flag)
            show_s_vals = [0 logspace(-5, -1, 9)]; % [0 logspace(-3, -1, 5)]; % New: plot only some values of s % logspace(-4, -1, 7)
        else
            show_s_vals = [0 logspace(-3, -1, 5)];
        end
        show_model_inds = intersect(model_inds, find(ismember(s_vec, show_s_vals)));
        switch plot_type{1}
            
            case 'coefficient_of_variation' % NEW! plot COV
                for j=[1 length(model_inds):-1:2]
                    coefficient_of_variation(j) = std(n_vec{model_inds(j)}) / mean(n_vec{model_inds(j)});
                end
                [cur_s_vec, sort_perm] = sort(s_vec(model_inds));
                semilogx(max(10^(-6), cur_s_vec), coefficient_of_variation(sort_perm), 'color', population_color_vec{i_d}, ...
                    'linewidth', 2); hold on;
                
                output_plot_file_name = '_coefficient_of_variation'; 
                
            case 'density' % ???
                for j=[1 length(model_inds):-1:2]
                    tmp_color_ind = mod_max(6-floor(mod(j, 10)/2), 5);
                    semilogx(n_hist_vec{model_inds(j)}, p_vec{model_inds(j)}, ...
                        [eric_color_vec(tmp_color_ind) my_symbol_vec{mod_max(j,2)}], 'linewidth', 2); hold on;
                end
                xlabel(str2title('f_null')); ylabel('Freq.');
                
                output_plot_file_name = '_num_alleles_dist';
                legend(legend_vec_only_s(model_inds([1 end:-1:2])));
                title(['Variation in f_{null}. ' model_str], 'fontsize', 14, 'fontweight', 'bold'); %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
            case 'inverse_cumulative_normalized'    % New: plot inverse cumulative normalized
                output_plot_file_name = '_inverse_cumulative_normalized';
                
                for j=[1 length(model_inds):-1:2]
                    if(ismember(model_inds(j), show_model_inds))
                        s_log = -log10(s_vec(model_inds(j)));
                        if(isinf(s_log))
                            s_log = 5.5;
                        end
                        jjj = j
                        tmp_color_ind = 6-floor(s_log); % mod_max(6-floor(mod(j, 10)/2), 5)
                        plot_x_vec = sort(n_vec{model_inds(j)}) ./ mean(n_vec{model_inds(j)}); % Plot normalized variation
                        plot_p_vec = (num_iters(i):-1:1) ./ num_iters(i);  cum_str = ' (inverse cumulative)'; % get cumulatives normalized
                        plot(plot_x_vec, plot_p_vec, 'color', main_selection_color_vec{tmp_color_ind}, ...
                            'linestyle', my_symbol_vec{mod_max(j,2)}, 'linewidth', 2);             hold on;
                    end
                end
                %                plot(repmat(1, 101, 1), 0:0.01:1, 'k', 'linewidth', 2);
                if(~SI_flag)
                    xlabel('f_{null} / E[f_{null}]', 'fontsize', 14, 'fontweight', 'bold');
                    ylabel('Freq. (Inverse Cumulative)', 'fontsize', 14, 'fontweight', 'bold');
                else
                    if(ismember(i_d, [3 4]))
                        ylabel('Freq. (Inverse Cumulative)');
                    end
                    if(ismember(i_d, [5 6]))
                        xlabel('f_{null} / E[f_{null}]');
                    end
                end
                xlim([0 30]); ylim([0 0.1]); %  at what fold-increase to stop? 10? 20?
                if(~SI_flag)
                    
                    if(i_d == 3) % Expansion1
                        h_leg = legend(cellstr(legend_vec_only_s(show_model_inds([1 end:-1:2]))')); % , 'WestOutside');
                        set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
                    end
                end
                title(model_str, 'interpreter', 'latex'); % 14, 'fontweight', 'bold'); % 'Variation in f_{null}. ' %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
                %                text(0.82, -0.0025, '1');
                
            case 'CAF_cumulative_log'
                show_s_vals = [0 logspace(-5, -1, 9)]; % [0 logspace(-3, -1, 5)]; % New: plot only some values of s % logspace(-4, -1, 7)
                show_model_inds = intersect(model_inds, find(ismember(s_vec, show_s_vals)));
                
                output_plot_file_name = '_CAF_cumulative_log';
                
                for j=[1 length(model_inds):-1:2]
                    if(ismember(model_inds(j), show_model_inds))
                        s_log = -log10(s_vec(model_inds(j)));
                        if(isinf(s_log))
                            s_log = 5.5;
                        end
                        jjj = j
                        tmp_color_ind = 6-floor(s_log); % mod_max(6-floor(mod(j, 10)/2), 5)
                        plot_x_vec = sort(n_vec{model_inds(j)}); % NOT normalized  ./ mean(n_vec{model_inds(j)}); % Plot normalized variation
                        plot_p_vec = (1:num_iters(i)) ./ num_iters(i);  cum_str = ' (cumulative)'; % get cumulatives normalized
                        semilogx(plot_x_vec, plot_p_vec, 'color', selection_color_vec{tmp_color_ind}, ...
                            'linestyle', my_symbol_vec{mod_max(j,2)}, 'linewidth', 2); hold on;
                    end
                end
                %                plot(repmat(1, 101, 1), 0:0.01:1, 'k', 'linewidth', 2);
                xlim([10^(-6.05) 1.1]); ylim([0 1]); %  at what fold-increase to stop? 10? 20?
                if(ismember(i_d, [5 6]))
                    xlabel('Combined Allele Frequency, $f_{null}$', 'interpreter', 'latex'); % , 'fontsize', 14, 'fontweight', 'bold');
                end
                if(ismember(i_d, [3 4]))
                    ylabel('Proportion of Genes with CAF below f_{null}'); % , 'interpreter', 'latex'); % , 'fontsize', 14, 'fontweight', 'bold');
                end
                title(model_str, 'interpreter', 'latex'); % , 'fontsize', 14, 'fontweight', 'bold'); % 'Variation in f_{null}. ' %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
                %                text(0.82, -0.0025, '1');
                
                
            case 'median_allele_frequency'
                output_plot_file_name = '_median_IAF_cumulative';
                for j=[1 length(model_inds):-1:2]
                    %                    if(ismember(model_inds(j), show_model_inds))
                    tmp_color_ind = mod_max(6-floor(mod(j, 10)/2), 5);
                    plot_x_vec = sort(median_vec{model_inds(j)}); % Plot normalized variation
                    plot_p_vec = (1:num_iters(i)) ./ num_iters(i);  cum_str = ' (median IAF)'; % get cumulatives normalized
                    semilogx(plot_x_vec, plot_p_vec, [eric_color_vec(tmp_color_ind) my_symbol_vec{mod_max(j,2)}], 'linewidth', 2);             hold on;
                    %                    end
                end
                xlabel('Median Allele Frequency', 'fontsize', 14, 'fontweight', 'bold');
                ylabel('Freq. (Inverse Cumulative)', 'fontsize', 14, 'fontweight', 'bold');
                
                % No legend (one legend outside)               legend(legend_vec_only_s(model_inds([1 end:-1:2])));
                title(['Variation in Median IAF. ' model_str], 'fontsize', 14, 'fontweight', 'bold'); %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
                
                
                %             case 's_estimator_variance' % NEW! here we plot the variance in the estimated s. Where we estimate: \hat{s} = mu / f
                %                 plot_x_vec = s_vec(model_inds); % One plot as function of s
                %                 plot_y_vec = std_vec;
                %
                %                 errorbar(plot_x_vec, plot_y_vec, std_vec);  % NEED TO MODIFY HERE !!! (NOT NEWEST VERSION!)
                
            case 's_confidence_interval' % show variation in estimateds: \hat{s} = \mu / f_C
                output_plot_file_name = '_s_confidence_interval';
                mu = 10^(-5); % rate used for simulations per gene
                
                plot_y_vec = s_vec(model_inds); plot_x_vec = zeros(length(model_inds), 1); plot_x_minus_std_vec = plot_x_vec;
                plot_x_plus_std_vec = plot_x_vec;
                s_ctr=1;
                for j=vec2row(model_inds)
                    plot_x_vec(s_ctr) = median(mu ./ n_vec{j});
                    plot_x_minus_std_vec(s_ctr) = quantile(mu ./ n_vec{j}, 0.025);
                    plot_x_plus_std_vec(s_ctr) = quantile(mu ./ n_vec{j}, 0.975);
                    s_ctr=s_ctr+1;
                end
                plot_x_plus_std_vec(isnan(plot_x_plus_std_vec)) = 1; %if no mutations observed, we estimate s=1
                loglog(plot_x_vec([1 end:-1:2]), plot_y_vec([1 end:-1:2]), 'linewidth', 2); hold on;
                loglog(plot_y_vec([1 end:-1:2]), plot_y_vec([1 end:-1:2]), '--', 'linewidth', 2); hold on;
                loglog(plot_x_minus_std_vec([1 end:-1:2]), plot_y_vec([1 end:-1:2]), 'k--', 'linewidth', 2);
                loglog(plot_x_plus_std_vec([1 end:-1:2]), plot_y_vec([1 end:-1:2]), 'k--', 'linewidth', 2);
                ylim([10^(-4) 10^(-1)]); xlim([min(plot_x_minus_std_vec) * 0.9 max(plot_x_plus_std_vec) * 1.1]); % set limits to fig
                for chosen_s = [10^(-2) 10^(-3)]
                    [~, chosen_s_plus_ind] = min(abs(plot_x_plus_std_vec - chosen_s)); % Plot example of confidence interval for one particular s
                    [~, chosen_s_minus_ind] = min(abs(plot_x_minus_std_vec - chosen_s)); % Plot example of confidence interval for one particular s
                    line([chosen_s chosen_s], [plot_y_vec(chosen_s_plus_ind) plot_y_vec(chosen_s_minus_ind)], 'color', 'r', 'linewidth', 2);
                end
                %                set(gca, 'xscale', 'log'); set(gca, 'yscale', 'log');
                xlabel('Estimated selection coefficient $\hat{s}$', 'fontsize', 14, 'fontweight', 'bold', 'interpreter', 'latex');
                ylabel('True selection coefficient $s$', 'fontsize', 14, 'fontweight', 'bold', 'interpreter', 'latex');
                title(['Distribution of Estimated s, $\hat{s}=\frac{\mu}{f_C}$. ' model_str], 'fontsize', 14, 'fontweight', 'bold', 'interpreter', 'latex'); %   Mean=' num2str(mean(f_null)) '.St.d.=' num2str(std(f_null))]);
                h_leg = legend({'Median $\hat{s}$', 'True $s$', '$95\%$ confidence'}, 4, 'interpreter', 'latex');
                set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
            case 'alpha_confidence_interval' % NEW (July 31st) - get also variance of estimator of alpha (using missense)
                
                
        end % switch plot type
        %        legend('boxoff'); % legned_vec_no_pop
        if((i_d == 6) || (~strcmp(plot_type{1}, 'coefficient_of_variation')))
            add_faint_grid(0.5);
        end
        
        switch plot_type{1}
            case {'CAF_cumulative_log', 'coefficient_of_variation'}
                one_plot_flag = 1;
            case 'inverse_cumulative_normalized'
                if(SI_flag)
                    one_plot_flag = 1;
                else
                    one_plot_flag = 0;
                end
            otherwise
                one_plot_flag = 0;
        end
        
        if(one_plot_flag)
            if(i_d == 6) % Expansion1
                
                switch plot_type{1}
                    case 'coefficient_of_variation'
                        h_leg = legend(pop_legend_vec, 2); % 'eastoutside'); % legend('boxoff'); % just legend
                        set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
                    otherwise                        
                        h_leg = legend(cellstr(legend_vec_only_s(show_model_inds([1 end:-1:2]))')); % , 'WestOutside');
                        pos_leg = get(h_leg, 'position'); set(h_leg, 'position', [pos_leg(1)+0.15 pos_leg(2)+0.21 pos_leg(3) pos_leg(4)]); % 275, -015
                        %                    set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
                        legend('boxoff'); % remove completely on side
                end
                
                orient landscape;
                my_saveas(gcf, fullfile(figs_dir, 'different_models/variation', ...
                    ['all_six_populations'  output_plot_file_name] ), {'epsc', 'pdf', 'jpg'});
            end
            
        else
            my_saveas(gcf, fullfile(figs_dir, 'different_models/ALL_variation', ...
                [unique_pop_vec{i_d}  output_plot_file_name] ), {'epsc', 'pdf', 'jpg'});
        end
        
    end % loop on plot type
end % loop on models




for i_s = 1:length(unique_s_vec) % Loop on s and plot separately for each s and all populations
    if(~ismember(unique_s_vec(i_s), [0 10.^[-1:-0.5:-5]]))
        continue;
    end
    for plot_type = {'inverse_cumulative_normalized'} % {'density', 'inverse_cumulative', 'inverse_cumulative_normalized'}
        s_inds = find(unique_s_vec(i_s) == s_vec);
        figure;
        bad_j = []; model_ctr = 1;
        for j=1:length(s_inds) % loop on models
            switch pop_vec{s_inds(j)}
                case {'expan1'} % , '2phase'} % 'expan2'
                    cur_symbol_vec = '-';
                case {'europ', 'ice3', 'finn3'} % , 'finn2'
                    cur_symbol_vec = '-';
                otherwise % case {'expan2', 'finn2'} % we don't do these
                    bad_j = [bad_j j];
                    continue;
            end
            plot_model_ctr = strmatch(pop_vec{s_inds(j)}, six_populations_vec)
            switch plot_type{1}
                case 'density'
                    plot_x_vec = n_hist_vec{s_inds(j)};
                    plot_p_vec = p_vec{s_inds(j)}; cum_str = '';
                    semilogx(plot_x_vec, plot_p_vec, [color_vec(model_ctr) cur_symbol_vec], 'linewidth', 2);             hold on;
                    xlabel(str2title('f_{null}'));
                case 'inverse_cumulative'
                    plot_x_vec = sort(n_vec{s_inds(j)});
                    plot_p_vec = (num_iters(i):-1:1) ./ num_iters(i);  cum_str = ' (inverse cumulative)'; % get inverse cumulatives
                    semilogx(plot_x_vec, plot_p_vec, [color_vec(model_ctr) cur_symbol_vec], 'linewidth', 2);             hold on;
                    xlabel(str2title('f_{null}'));
                case 'inverse_cumulative_normalized'
                    plot_x_vec = sort(n_vec{s_inds(j)}) ./ mean(n_vec{s_inds(j)}); % Plot normalized
                    plot_p_vec = (num_iters(i):-1:1) ./ num_iters(i);  cum_str = ' (inverse cumulative)'; % get cumulatives normalized
                    plot(plot_x_vec, plot_p_vec,  cur_symbol_vec, ...
                        'color', population_color_vec{plot_model_ctr}, 'linewidth', 2); hold on;
                    xlabel('f_{null} / E[f_{null}]', 'fontsize', 14, 'fontweight', 'bold');                    
            end
            model_ctr = model_ctr+1;
        end
        ylabel(['Freq.' cum_str], 'fontsize', 15, 'fontweight', 'bold');
        title(['Variation in f_{null}. s=' s_vec_str{s_inds(1)}], 'fontsize', 15, 'fontweight', 'bold'); %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
        h_leg = legend(legend_vec_no_s(s_inds(setdiff(1:length(s_inds), bad_j))), 'fontsize', 12, 'fontweight', 'bold');
        set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
        set(gca, 'fontsize', 15);
        switch plot_type{1}
            case {'density', 'inverse_cumulative'}
                y_lim = ylim(gca); ylim([y_lim(1) min(y_lim(2), 10000)]);
            case 'inverse_cumulative_normalized' % This is the one we use!
                xlim([0.0 30]);
                ylim([0 0.1]); % Set minimal 'interesting' frequency
                %                x_lim = xlim(gca); xlim([0.99 x_lim(2)]);
                x_tick = get(gca, 'XTick');
                if(x_tick(1) > 1)
                    x_tick = [1 x_tick];
                    set(gca, 'XTick', x_tick);
                end
        end
        add_faint_grid(0.5);
        my_saveas(gcf, fullfile(figs_dir, plot_type{1}, ...
            ['num_alleles_' plot_type{1} '_mu_1_7_10_6_s_' strrep(strrep(strrep(s_vec_str{s_inds(1)}, '^', '_'), '{', ''), '}', '_')] ), {'epsc', 'pdf'});
    end % loop on plot type
    % close all;
    
    % New: also make boxplot
    figure; boxplot(0.2.*cell2mat(n_vec(s_inds)')', 'labels', pop_vec(s_inds));     set(gca,'yscale','log');
    hold on;
    plot(mean(0.2.*cell2mat(n_vec(s_inds)')'), '*k', 'markersize', 10);
    title(['Variation in n-null. s=' s_vec_str{s_inds(1)}]);
    my_saveas(gcf, fullfile(figs_dir, 'boxplot', ...
        ['num_alleles_boxplot_s_' strrep(strrep(strrep(s_vec_str{s_inds(1)}, '^', '_'), '{', ''), '}', '_')] ), {'epsc', 'pdf'});
    
    
    
    
    
    % New: plot quantiles as 'my boxplots'
    figure;
    quantiles_vec = [0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99];
    quantile_mat = zeros(length(quantiles_vec), length(s_inds));
    for j=1:length(quantiles_vec)
        quantile_mat(j,:) = my_quantile(0.2.*cell2mat(n_vec(s_inds)')', quantiles_vec(j));
    end
    median_vec = quantile_mat(5,:);
    set(gca,'yscale','log'); hold on;
    
    plot(mean(0.2.*cell2mat(n_vec(s_inds)')'), '*k', 'markersize', 12);
    quantile_color_vec = 'bgmrk';
    for j=5:-1:1 % loop on quantiles length(quantiles_vec)/2 % plot boxes
        if(j==5)
            line_width = 3;
        else
            line_width = 1.25;
        end
        for k=1:length(s_inds) % loop on different models
            rectangle('position', [k-0.25, max(0.000001, quantile_mat(j, k)), ...
                0.5, max(0.0000000001, quantile_mat(10-j, k) - quantile_mat(j, k))], 'edgecolor', quantile_color_vec(j), 'linewidth', line_width); hold on;
        end
        plot(1, quantile_mat(j, 1), 'color', quantile_color_vec(j), 'linewidth', line_width); hold on; % for legend
    end
    plot(mean(0.2.*cell2mat(n_vec(s_inds)')'), '*k', 'markersize', 12);
    set(gca, 'XTick', 1:length(s_inds))
    set(gca, 'XTicklabel', pop_vec(s_inds), 'fontsize', 14);
    legend({'mean', 'median', '(25%,75%)', '(10%,90%)', '(5%,95%)', '(1%,99%)'}, 'location',  'South');    legend('boxoff')
    title(['Variation in n-null. s=' s_vec_str{s_inds(1)}], 'fontsize', 16);
    ylim([10^(-6) 10^(-1)]);
    my_saveas(gcf, fullfile(figs_dir, 'boxplot', ...
        ['num_alleles_my_boxplot_s_' strrep(strrep(strrep(s_vec_str{s_inds(1)}, '^', '_'), '{', ''), '}', '_')] ), {'epsc', 'pdf'});
    
    
    % % %     figure; % New: plot quantiles in a different way: cumulative curves for each model
    % % %     for k=1:length(s_inds) % loop on different models
    % % %         cur_model_f_vec = 0.2.*cell2mat(n_vec(s_inds)')';
    % % %         hold on; plot(sort(cur_model_f_vec), (1:length(cur_model_f_vec)) ./ length(cur_model_f_vec), color_vec(j));
    % % %     end
    % % %     title(['Variation in f_{null}. s=' s_vec_str{s_inds(1)}], 'fontsize', 16);
    % % %     xlabel('f_{null}'); ylabel('Cumulative frequency');
    % % %     my_saveas(gcf, fullfile(figs_dir, 'boxplot', ...
    % % %         ['num_alleles_cumulative_s_' strrep(strrep(strrep(s_vec_str{s_inds(1)}, '^', '_'), '{', ''), '}', '_')] ), {'epsc', 'pdf'});
    
    
end % loop on s


% New: make


% New: make quantile plots for each model, for different values of s
for pop_groups = { {'equil', 'expan1', 'expan2', '2phase'}, {'europ', 'ice', 'finn1', 'finn2'} } % , {'varsel1', 'varsel2'}, ...
    %%%        {'finn1_bneck', 'finn1_nobneck', 'finn2_bneck', 'finn2_nobneck'}, {'ice_bneck', 'ice_nobneck'} }; % new! add finland and iceland
    if(ismember('equil', pop_groups{1}))
        pop_str = 'expansion';
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
    end
    figure;
    model_ctr=1;
    for i_d = 1:num_models % Loop on models and plot separately for each population
        model_inds = strmatch(unique_pop_vec{i_d}, pop_vec);
        model_str = get_nice_population_names(unique_pop_vec{i_d});
        
        if(ismember(unique_pop_vec{i_d}, pop_groups{1}))
            subplot(2,2,model_ctr);
            
            %             quantiles_legend_vec = plot_quantiles(cell2mat(n_hist_vec(model_inds )'), cell2mat(p_vec(model_inds )'), ...
            %                 max(10^(-6), s_vec(model_inds)), 1);
            
            quantiles_legend_vec = plot_quantiles([], 0.2.*cell2mat(n_vec(model_inds)'), ... % divide by 5 to get #LOF alleles.
                max(10^(-6), s_vec(model_inds)), 1);
            
            xlabel(str2title('Selection Coefficient s')); ylabel('Number of LOF alleles per chromosome');
            title(str2title(['Variation in f_{null}. ' model_str])); %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
            ylim([10^(-4) 10^1]);
            add_faint_grid(0.7);
            model_ctr=model_ctr+1;
        end % found model
    end % loop on models
    
    h_l = legend({quantiles_legend_vec}, 3, 'fontsize', 10); % 'location', 'east');
    pos_l = get(h_l, 'position'); set(h_l, 'position', [pos_l(1)+0.275 pos_l(2)-0.15 pos_l(3) pos_l(4)]); % change legend position
    legend('boxoff'); % +plot_type*(~figs_for_paper_flag));
    
    orient landscape;
    my_saveas(gcf, fullfile(figs_dir, 'different_models', [pop_str  '_num_alleles_quantiles'] ), {'epsc', 'pdf'});
end % loop on pop-groups


XXX = 234234;













