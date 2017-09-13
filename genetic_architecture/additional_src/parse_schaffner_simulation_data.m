% function parse_schaffner_simulation_data() % read output of Schaffner simulations

%demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/sim_results_set2';
%demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/files_var_eric'; % New Version !!
demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/EuropeFixed/files_europe_fixed/files_var_eric'; % New Version !! Europe fixed

figs_dir = '../../common_disease_model/figs/EyreWalker/new_eric/all_models/empirical_distribution';
allele_freq_file = fullfile(demographic_sims_dir, 'fbin.txt');

parse_mean = 0; % Look at summary files
parse_all_distribution = 1; % Look at individual simulations files
parse_age_data = 0; parse_age_data_to_mat=0; % NEW! parse allelic age files
plot_all_age_simulations = 0;

if(parse_mean)
    sims_files = union( GetFileNames(fullfile(demographic_sims_dir, '*.all')), ...
        GetFileNames(fullfile(demographic_sims_dir, '*.old')));
    sims_files = union(sims_files, GetFileNames(fullfile(demographic_sims_dir, '*.new')));
    sims_files = union(sims_files, GetFileNames(fullfile(demographic_sims_dir, '*.*bneck')));
    
    x_vec = ReadDataFile(allele_freq_file); x_vec = x_vec.derived_freq; % Note: allele freq. is in logarithmic coordinates.
    s_vec =  [10.^(-1:-0.5:-5) 0]; % , s=10^-1.5, s=10^-2, ... s=10^-5, s=0); % hard-coded selection based on Steve's mail. (discrapancy: 11 or 10 s?)
    
    my_mkdir(fullfile(demographic_sims_dir, 'mat'));
    for i=1:length(sims_files)
        cur_output_file =  [fullfile(demographic_sims_dir, 'mat', strrep(sims_files{i}, '.', '_')) '.mat'];
        p_vec = load(fullfile(demographic_sims_dir, sims_files{i}));
        if(strmatch('cumul', sims_files{i}))
            num_alleles_per_chrom = p_vec(:,end);
            save(cur_output_file, 'x_vec', 'p_vec', 's_vec', 'num_alleles_per_chrom');
        else
            num_sims = sum(p_vec,2); % This is num. of simulations???
            save(cur_output_file, 'x_vec', 'p_vec', 's_vec', 'num_sims');
        end
    end % loop on files
end % if parse mean


demographic_sims_all_distribution_dirs = {'../../common_disease_model/data/schaffner_simulations/out_var', ...
    '../../common_disease_model/data/schaffner_simulations/out_var2', ...
    '~sfs/working/rare_eric/out_var', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/files_europe_fixed/out_var', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/fine_s_var/out_var', ...  % NEW! fine-s simulation
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/fine_s_var_f_disruptive_1_7_10_6/out_var', ... % New Version! with mutation rate of 1.7*10^(-6) and not 10^(-5) per gene (only disruptive alleles) !
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/correct_fine_s_high_res', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/three_mu_fine_s/mu_10_5', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/three_mu_fine_s/mu_5_10_5', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/three_mu_fine_s/mu_1_7_10_6', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/stop_at_bottleneck_fine_s'};
%        '/fg/hapmaphet/out_var'}; % Added Europe  in broad unix
demographic_sims_all_distribution_output_mat_dirs = {'../../common_disease_model/data/schaffner_simulations/out_var/mat', ...
    '../../common_disease_model/data/schaffner_simulations/out_var2/mat', ...
    '../../common_disease_model/data/schaffner_simulations/out_var3/mat', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/files_europe_fixed/out_var/mat', ... % Added Europe
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/fine_s_var/out_var/mat', ... % NEW! fine-s simulations
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/fine_s_var_f_disruptive_1_7_10_6/out_var/mat', ... % New Version! with mutation rate of 1.7*10^(-6) and not 10^(-5) per gene (only disruptive alleles) !
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/correct_fine_s_high_res/mat', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/three_mu_fine_s/mu_10_5/mat', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/three_mu_fine_s/mu_5_10_5/mat', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/three_mu_fine_s/mu_1_7_10_6/mat', ...
    '../../common_disease_model/data/schaffner_simulations/EuropeFixed/stop_at_bottleneck_fine_s/mat'};
    

if(parse_all_distribution) % NEW! parse the entire distribution (not just means)
    preprocess_flag = [0 0 0 0 0 0 0 0 1 0 0]; % run only new fine variation for s. 5 is mu=10^{-5} and 7 is for mu=1.7*10^{-6}
    
    for i_dir = 9:9 % 5:5 % 8:8 % 7:7 % 5:5 need to copy new data !!! 5:length(demographic_sims_all_distribution_dirs) % loop on different runs
        i_dir_is = i_dir
        dist_file_names = GetFileNames(fullfile(demographic_sims_all_distribution_dirs{i_dir}, '*.all'), 1)
        num_files_to_preprocess = length(dist_file_names)
        
        if(preprocess_flag(i_dir)) % preprocess files - convert to .mat
            run_preprocess = 999
            for i=1:length(dist_file_names)
                sprintf('Run preprocess-file %ld out of %ld file: %s', i, length(dist_file_names), dist_file_names{i})
                run_dir = i_dir
                s_distribution = load( dist_file_names{i}); % load and convert to .mat file
                save(fullfile(demographic_sims_all_distribution_output_mat_dirs{i_dir}, ...
                    remove_dir_from_file_name(file_name_to_mat(dist_file_names{i}))), 's_distribution'); %
            end
        end
        plot_two_class_empirical_distributions(demographic_sims_all_distribution_output_mat_dirs{i_dir}, ...
            allele_freq_file, fullfile(figs_dir, ['var' num2str(i_dir)])); % parse and plot empirical distributions
        
    end % loop on dirs
end



if(parse_age_data) % New: Parse age data
    for i_dir = 7:7
        
        if(plot_all_age_simulations) % plot results
            orange = [1 0.6 0.1]; % set color
            population_color_vec = {'k', orange, 'r', 'g', 'b', 'c'};
            eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
            my_symbol_vec = {'--', '-'}; % flip ordering (set integer powers as solid lines)
            selection_color_vec = {'k', 'b', 'g', orange, 'r'}; % replace yellow with orange
            
            figure;
            load(fullfile(demographic_sims_all_distribution_dirs{i_dir}, 'mat', 'all_age_simulations.mat')); % load simulations data
            populations_vec = {'equil', 'europ', 'expan1', 'finn1', '2phase', 'ice'};
            show_s_vec = [1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 0.0];
            [~, selection_inds] = intersect_all(s_vec, show_s_vec);
            
            
            for i_pop = 1:length(populations_vec) % loop on populations
                show_s_vec = [1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 0.0];
                pop_inds = strmatch(populations_vec{i_pop}, pop_vec);
                pop_selection_inds = intersect(pop_inds, selection_inds);
                if(~ismember(1, s_vec(pop_selection_inds)))  % hack (missing steven's simulations
                    pop_selection_inds = [pop_selection_inds(1) ...
                        intersect(pop_inds, find(s_vec == 1.1)) pop_selection_inds(2:end)]; % add 1.1
                    show_s_vec(1)=1.1;
                end
                if(~ismember(2.5, s_vec(pop_selection_inds)))  % hack (missing steven's simulations
                    pop_selection_inds = [pop_selection_inds(1:4) ...
                        intersect(pop_inds, find(s_vec == 2.4)) pop_selection_inds(5:end)]; % add 1.1
                    show_s_vec(4)=2.4;
                end
                subplot(3,2,i_pop);
                for j=length(show_s_vec):-1:1
                    [~, cur_plot_ind] = intersect(s_vec(pop_selection_inds), show_s_vec(j));
                    cur_plot_ind = pop_selection_inds(cur_plot_ind);
                    if(isempty(cur_plot_ind))
                        continue;
                    end
                    
                    allele_age_cum{cur_plot_ind} = cumsum(allele_age_hist{cur_plot_ind}); % compute NORMALIZED cumulative
                    allele_age_cum{cur_plot_ind} = allele_age_cum{cur_plot_ind} ./ allele_age_cum{cur_plot_ind}(end);
                    semilogx(allele_age_bins{cur_plot_ind}, allele_age_cum{cur_plot_ind}, ...
                        'color', selection_color_vec{ceil((11-j)/2)}, 'linestyle', my_symbol_vec{mod_max(j,2)}, 'linewidth', 2);
                    hold on;
                end % loop on selection coefficient
                xlim([0 10^5]);
                if(i_pop >= 5)
                    xlabel('Age (num. generations, k)');
                end
                if((i_pop == 3) || (i_pop == 4))
                    ylabel('Proportion of alleles born in last k generations');
                end
                %                set(gca, 'fontsize', 14);
                title(get_nice_population_names(populations_vec{i_pop}));
                if(i_pop == 6) % legend
                    legend_vec = cellstr([repmat('s=10^{-', 10, 1) num2str(show_s_vec(end:-1:1)') repmat('}', 10, 1)]);
                    legend_vec = strrep_cell(legend_vec, ' ', '');
                    legend_vec = strrep_cell(legend_vec, '10^{-0}', '0');
                    h_leg = legend(legend_vec, 3);
                    pos_l = get(h_leg, 'position'); set(h_leg, 'position', [pos_l(1)+0.2835 pos_l(2)-0.1435 pos_l(3) pos_l(4)]);
                    set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
                end
                
                add_faint_grid(0.5);
            end % loop on populations
            orient landscape;
            
            my_saveas(gcf, fullfile(figs_dir, 'age_cumulative_populations'), ...
                {'epsc', 'pdf', 'jpg'});
            
            
            % NEW: Plot median allele frequency for each allele
            close all;
            num_files = length(allele_age_hist);
            for i=1:num_files
                age_median_vec(i) = median_hist(allele_age_bins{i}, allele_age_hist{i});
            end
            fig_var = [];
            for i_fig=1:1 % 2
                if(i_fig == 1) % First is main figure 
                    fig_var_X = figure(i_fig);
                else % second is inset figure 
                    fig_var_Y = figure(i_fig);
                end
                i_pop_leg=1;
                for i_pop = [1 3 5 2 4 6] % 1:length(populations_vec) % loop on populations
                    pop_inds = strmatch(populations_vec{i_pop}, pop_vec);
                    tmp_s_vec = 10.^(-s_vec(pop_inds)); tmp_s_vec(tmp_s_vec==1)=0;
                    [~, sort_perm] = sort(tmp_s_vec);
                    semilogx(tmp_s_vec(sort_perm), age_median_vec(pop_inds(sort_perm)), 'color', population_color_vec{i_pop_leg}, 'linewidth', 2);
                    hold on;
                    i_pop_leg=i_pop_leg+1;
                end
                if(i_fig==1)
                    xlabel('Selection coefficient, s', 'fontsize', 15); ylabel('Median age (num. generations)', 'fontsize', 15);
                else
                    xlim([10^(-3) 10^(-1)]);
                end
                set(gca, 'fontsize', 14);
                y_ticks = get(gca, 'ytick');
                y_tick_labels = num2str(y_ticks');
                set(gca, 'yticklabel', y_tick_labels);
                y_tick_labels = get(gca, 'yticklabel')
                
                
                h_leg = legend(get_nice_population_names(populations_vec([1 3 5 2 4 6])), 1, 'fontsize', 14);
                set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
                
                add_faint_grid(0.5);
                orient landscape;
                my_saveas(gcf, fullfile(figs_dir, 'age_median_populations'), ...
                    {'epsc', 'pdf', 'jpg'});
            end % loop on figure;

            
%            [h_m h_i]=inset(fig_var_X,fig_var_Y, 0.72); % skip inset for             now 
%            h_leg = legend(h_i, get_nice_population_names(populations_vec([1 3 5 2 4 6])), 1, 'fontsize', 14);
%            set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
            
            add_faint_grid(0.5);
            orient landscape;
            
            my_saveas(gcf, fullfile(figs_dir, 'age_median_populations_with_inset'), ...
                {'epsc', 'pdf', 'jpg'});
            
            
        else % here parse age data
            
            look_in_dir = demographic_sims_all_distribution_dirs{i_dir}
            age_file_names = GetFileNames(fullfile(demographic_sims_all_distribution_dirs{i_dir}, '*.age'), 1)
            num_files = length(age_file_names);
            allele_age_hist = cell(num_files,1); allele_age_bins = cell(num_files,1);
            s_vec = zeros(num_files,1); pop_vec = cell(num_files,1);
            age_file_names_mat = cell(num_files, 1);
            if(~parse_age_data_to_mat)
                age_file_names_mat = GetFileNames(fullfile(demographic_sims_all_distribution_dirs{i_dir}, 'mat', 'age', '*.age*mat'), 1)
                num_files = length(age_file_names_mat);
            end
            run_preprocess_age_data = 999
            for i=1:num_files
                run_file = i
                if(~parse_age_data_to_mat)
                    age_file_names{i} = remove_suffix_from_file_name(strrep(age_file_names_mat{i}, '\mat', ''));
                end
                
                parse_job_str = ['parse_age_data_internal(''' age_file_names{i} ''')'];
                tmp_file_name = remove_dir_from_file_name(age_file_names{i});
                pop_vec{i} = str2word('.', tmp_file_name, 1);
                dot_ind = strfind(tmp_file_name, '.');
                s_vec(i) = str2nums(tmp_file_name(dot_ind(1)+1:end));
                
                if(parse_age_data_to_mat)
                    age_file_names_mat{i} = fullfile(demographic_sims_all_distribution_dirs{i_dir}, 'mat', 'age', ...
                        [remove_dir_from_file_name(age_file_names{i}) '.mat']);% Create histograms
                    if(~exist(age_file_names_mat{i}, 'file'))
                        SubmitMatlabJobToFarm(parse_job_str, ...
                            fullfile(demographic_sims_all_distribution_dirs{i_dir}, 'out', ...
                            [remove_dir_from_file_name(age_file_names{i}) '.out']), 'hour', ...
                            [], [], 8); % NEW! add memory for new big files!!!
                    end
                else % here already in .mat file
                    load(age_file_names_mat{i});
                    [allele_age_hist{i} allele_age_bins{i}] = ...
                        weighted_hist(cell2vec(allele_age_vec), cell2vec(allele_freq_vec), ...
                        1:max(max_cell(allele_age_vec(~isempty_cell(allele_age_vec)))), 1); % get histogram
                end % if to-mat
            end % loop on age files
            if(~parse_age_data_to_mat)
                save(fullfile(demographic_sims_all_distribution_dirs{i_dir}, 'mat', 'all_age_simulations.mat'), ...
                    'allele_age_hist', 'allele_age_bins', 'age_file_names', 'age_file_names_mat', 's_vec', 'pop_vec');
            end
            
            
        end % plot all simulations results or parse data
        
    end % loop on directory
end % if parse age data


