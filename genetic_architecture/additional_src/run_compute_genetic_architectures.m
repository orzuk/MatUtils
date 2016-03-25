% Simulate many architectures and check if they give desired statistics
%
% Input:
% architecture_str_vec - string representing architectures type to run
% action_mode - run computations (default) or just collect results
% in_matlab_flag - run inside matlab (default) or submit jobs to farm
% N_vec - vector of #loci
% f_vec - vector of minor-allele-frequencies
% disease_type_str - string representing the disease
% queue_str - where to submit job to
%
% Output:
% good_architectures - structure with output architecture
% good_architectures_file - file where architectures were saved
% good_architectures_plot_file - file where figures were saved
% good_architectures_latex_file - file where summary statistics were saved in latex format
%
function [good_architectures good_architectures_file good_architectures_plot_file good_architectures_latex_file] = ...
    run_compute_genetic_architectures(architecture_str_vec, action_mode, in_matlab_flag, ...
    N_vec, f_vec, disease_type_str, queue_str)

intervals_struct = disease_type_to_constraint(disease_type_str);

if(~exist('action_mode', 'var') || isempty(action_mode))
    action_mode = 1;
end
if(~exist('in_matlab_flag', 'var') || isempty(in_matlab_flag))
    in_matlab_flag = 1;
end
AssignGeneralConstants;
AssignStatsConstants;
[machine machine_delim html_outdir] = get_machine_type();
working_location = get_working_location();
if(~exist('queue_str', 'var') || isempty(queue_str))
    queue_str = 'broad'; % 'broad' some archs. need lots of memory (?)
end
compute_method_flag_vec = {'analytic', 'enumerate'}; %compute_method_flag = 'analytic'; % 'enumerate'; % 'analytic'; % 'enumerate'; % 'analytic'; % 'enumerate';            ttt_sim = cputime;
num_rep_per_genotype = [1 5];
dispose_prob_frac = 0; % ignore this part of the probability space
iters = 111; % number of functions to test for each type of architecture
plot_flag = 0; % if to plot figures during work
power_flag = 1; % if to compute power to detect each architecture in a GWAS
params_struct = ...
    set_architecture_params(max(N_vec), 'sigmoid', ...
    intervals_struct.h_interval, intervals_struct.freq_interval, intervals_struct.penetrance_interval, ...
    max(N_vec), 1, iters);

switch machine
    case UNIX
        good_architectures_dir = '/seq/orzuk/common_disease_model';
    case PC
        switch working_location
            case 'BROAD'
                good_architectures_dir = 'T:\common_disease_model';
            case 'HOME' % at home
                good_architectures_dir = 'c:/research/common_disease_model/';
        end
end

if(iscell(architecture_str_vec)) % input must be cell array for splitting!
    num_architectures = length(architecture_str_vec);
    if(action_mode) % run jobs
        good_architectures_file = []; good_architectures_plot_file = []; good_architectures_latex_file = []; % output dummy variables
        for i=1:num_architectures % loop on architectures
            for j = 1:length(N_vec) % loop on N
                if(in_matlab_flag)
                    run_computation_in_matlab = 999
                    good_architectures = run_compute_genetic_architectures(architecture_str_vec{i}, 1, 1, N_vec(j), f_vec, ...
                        disease_type_str, queue_str);
                else
                    good_architectures = [];
                    SubmitMatlabJobToFarm(['run_compute_genetic_architectures(''' architecture_str_vec{i} ...
                        ''', 1, 1, ' num2str(N_vec(j)) ', [ ' num2str(f_vec) '], ''' ...
                        disease_type_str ''', ''' queue_str ''');'], ...
                        fullfile(good_architectures_dir, ['out/run_architectures_' architecture_str_vec{i} ...
                        '_N_' num2str(N_vec(j)) '.out']), ...
                        queue_str);
                end
            end
        end
        return; % no saving of one file at the end (??)
    else % collect results
        good_architectures = []; plot_flag = 1; % plot all (not ONLY global plots (not architecture specific))
        for i=1:num_architectures % load results files and concatenate them
            for j = 1:length(N_vec) % loop on N
                good_architectures_file = ...
                    get_architecture_file_names(good_architectures_dir, architecture_str_vec{i}, N_vec(j));
                if(exist(good_architectures_file, 'file')) % some might be empty (no architectures)
                    OLD = load(good_architectures_file);
                    for k=1:length(OLD.good_architectures)
                        OLD.good_architectures(k).index = k;
                        if(~isfield(OLD.good_architectures(k), 'plot_name_arch'))
                            OLD.good_architectures(k).plot_name_arch = ...
                                arch_name_to_plot_name(OLD.good_architectures(k).arch);
                        end
                    end
                    good_architectures = [good_architectures OLD.good_architectures];
                end
            end
        end
        if(isempty(good_architectures))
            didnt_find_any_good_architecture = 99999
            return;
        end
        N=1;
        for i=1:length(good_architectures);
            N = max(N, good_architectures(i).N);
        end
        [good_architectures_file good_architectures_plot_file good_architectures_latex_file] = ...
            get_architecture_file_names(fullfile(good_architectures_dir, disease_type_str), ...
            [disease_type_str '_model_all_filtered_max'], N);
        valid_inds = ...
            filter_good_architectures(good_architectures, ...
            intervals_struct.h_interval, intervals_struct.h_add_interval, ...
            intervals_struct.freq_interval, intervals_struct.ratio_interval, ...
            intervals_struct.lods_interval, intervals_struct.penetrance_interval); % ...       [],     plot_flag);    % filter (some old ones may not be good for our interval)
        %        good_architectures = good_architectures(valid_inds);
        good_architectures = ...
            get_best_valid_architecture(good_architectures, valid_inds, ...
            params_struct, good_architectures(1).f_vec, good_architectures(1).arch, 1);
        
        %         if(plot_flag) % plot good architectures
        %             plot_architecture_statistics(good_architectures, ...
        %                 params_struct, good_architectures_plot_file);
        %         end
        %         save_architectures_in_latex(good_architectures, ...
        %             good_architectures_latex_file, ...
        %             good_architectures_file, good_architectures_plot_file); % save to nice latex tables format
    end % if action mode
else % here it's not cell - do the whole work
    % iters_is = iters
    % N_is = max(N_vec)
    % good_dir_is = good_architectures_dir
    k_in_clause_vec = [aliquotparts(N_vec / 2) * 2 N_vec]; % N must be even. Take all divisors
    num_freqs = length(f_vec);
    N_vec = repmat(N_vec, 1, num_freqs * length(k_in_clause_vec));
    f_vec = vec2row(mat_into_vec(repmat(f_vec, length(k_in_clause_vec), 1)));
    k_in_clause_vec = repmat(k_in_clause_vec, 1, num_freqs);
    num_clauses_vec = N_vec ./ k_in_clause_vec;
    num_models = length(N_vec);
    [good_architectures_file good_architectures_plot_file good_architectures_latex_file] = ...
        get_architecture_file_names(good_architectures_dir, architecture_str_vec, max(N_vec));
    
    if(exist(good_architectures_file, 'file'))
        OLD = load(good_architectures_file); % load old architectures
    end
    total_models_loop_time = cputime;
    interesting_vec = zeros(num_models,1);
    if(iscell(architecture_str_vec))
        architecture_str = architecture_str_vec{1}
    else
        architecture_str = architecture_str_vec
    end
    % good_architectures = struct([]); % good_architectures = cell(num_models,1);
    for k=1:num_models % loop on parameters
        do_model = k
        total_models = num_models
        N = N_vec(k); k_in_clause = k_in_clause_vec(k); num_clauses = num_clauses_vec(k);
        f = repmat(f_vec(k), 1, N);
        params_struct = ...
            set_architecture_params(N, architecture_str, ...
            intervals_struct.h_interval, intervals_struct.freq_interval, ...
            intervals_struct.penetrance_interval, ...
            k_in_clause, num_clauses, iters);
        compute_method_flag = compute_method_flag_vec{1};
        switch compute_method_flag
            case 'enumerate'
                [x_vec p_x_vec x_ind_vec x_ind_mat] = ...
                    initilize_x_vec_constants(N, dispose_prob_frac, f);
            otherwise
                x_vec = []; p_x_vec = []; x_ind_vec = []; x_ind_mat = [];
        end
        [candidate_architectures mu V v_marginal v_environment v_genetic V_helper ...
            v_marginal_explained v_additive_explained v_pairwise v_pairwise_explained ...
            h h_add h_liability h_add_fitted h_mult_fitted ...
            H_from_twins h_add_from_twins h_liability_from_twins h_liability_from_twins_ADE ...
            penetrance lods_ratio_marginal p_z_x_marginal, ...
            lods_ratio_pairwise, mu_pairwise, ...
            mu_given_k_ones, mu_given_k_ones_std, ...
            relative_risk, family_tree, family_risk, ...
            architecture_formula full_calc_vec] = ...
            compute_architecture_statistics_master(f, architecture_str, params_struct, ...
            compute_method_flag, num_rep_per_genotype(1), x_vec, p_x_vec, x_ind_mat, ...
            intervals_struct.h_interval, intervals_struct.h_add_interval, ...
            intervals_struct.ratio_interval, intervals_struct.freq_interval, ...
            intervals_struct.penetrance_interval, intervals_struct.lods_interval);
        candidate_architectures = struct('mu', mu, 'h', h, 'h_add', h_add, ...
            'h_liability', h_liability, ... % New! add the liability version of h
            'H_from_twins',  H_from_twins, ...
            'h_add_from_twins', h_add_from_twins, ...
            'h_liability_from_twins', h_liability_from_twins, ...
            'h_liability_from_twins_ADE', h_liability_from_twins_ADE, ...
            'V', V, 'v_marginal', v_marginal, 'v_environment', v_environment, ...
            'v_genetic', v_genetic, 'v_marginal_explained', v_marginal_explained, ...
            'v_additive_explained', v_additive_explained, 'v_pairwise', v_pairwise, ...
            'v_pairwise_explained', v_pairwise_explained, 'penetrance', penetrance, ...
            'lods_ratio_marginal', lods_ratio_marginal, ...
            'lods_ratio_pairwise', lods_ratio_pairwise, ...
            'p_z_x_marginal', p_z_x_marginal, 'mu_pairwise', mu_pairwise, ...
            'mu_given_k_ones', mu_given_k_ones, ...
            'mu_given_k_ones_std', mu_given_k_ones_std, ...
            'relative_risk', relative_risk, ...
            'family_tree', family_tree, ...
            'family_risk', family_risk, ...
            'architecture_formula', vec2column(architecture_formula), ...
            'full_calc_vec', full_calc_vec);
        [valid_inds valid_inds_mat interesting_vec(k)] = ...
            filter_good_architectures(candidate_architectures, ...
            intervals_struct.h_interval, intervals_struct.h_add_interval, ...
            intervals_struct.freq_interval, intervals_struct.ratio_interval, ...
            intervals_struct.lods_interval, intervals_struct.penetrance_interval); % ...       [],     plot_flag);
        good_architectures(k) = ...
            get_best_valid_architecture(candidate_architectures, valid_inds, ...
            params_struct, f, architecture_str, interesting_vec(k));
        if(k == 1) % get fields going to table only in the first time
            num_fields_in_table = good_architectures(k).num_fields_in_table;
        end
    end % loop on model parameters
    
    total_models_loop_time = cputime - total_models_loop_time
    if(~isempty(find(interesting_vec)))
        xxxx = 1241234124
    end
    good_architectures = good_architectures(find(interesting_vec)); % remove the ones not satisfying criteria
    
    if(power_flag) % compute power for each architecture (time consuming)
        total_power_calculation_time = cputime;
        for i=1:length(good_architectures)
            good_architectures(i).power_marginal  = ...
                compute_association_power(good_architectures(i).p_z_x_marginal, ...
                params_struct.n_samples_vec, [], params_struct.power_alpha, ...
                params_struct.power_iters*10, ... % we can afford more iters in marginal (faster)
                params_struct.power_test_type, params_struct.power_test_stat, ...
                'case-control');
            good_architectures(i).power_epistasis = ... % compute epistasis power: takes a long time ...
                compute_association_power(good_architectures(i).p_z_x_pairwise, ...
                params_struct.n_samples_vec, [], params_struct.power_pairwise_alpha, ...
                params_struct.power_iters, ...
                'epistasis', params_struct.power_pairwise_test_stat, ...
                'case-control');
        end
        total_power_calculation_time = cputime - total_power_calculation_time
    end
    
    if(exist(good_architectures_file, 'file')) % add old architectures (don't want to lose them)
        OLD = load(good_architectures_file);
        save_architectures_in_latex(good_architectures, ... % run this once to fill fields in new architectures list
            good_architectures_latex_file, ...
            good_architectures_file, good_architectures_plot_file); % save to nice latex tables format
        load(good_architectures_file);
        try
            good_architectures = [good_architectures OLD.good_architectures];
        catch exception % not possible to unite: take new one unless it's empty
            if(isempty(good_architectures))
                good_architectures = OLD.good_architectures;
            end
        end
    end % exists old architectures file 
end % if cell for architecture_str

save_architectures_in_latex(good_architectures, ...
    good_architectures_latex_file, ...
    good_architectures_file, good_architectures_plot_file, disease_type_str); % save to nice latex tables format again (with OLDs)
saved_to_latex = good_architectures_latex_file
if(plot_flag) % plot good architectures
    params_struct.plot_flag = plot_flag;
    params_struct.disease_type_str = disease_type_str;
    plot_architecture_statistics(good_architectures, ...
        params_struct, good_architectures_plot_file);
end


