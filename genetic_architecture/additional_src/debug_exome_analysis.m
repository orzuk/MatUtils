%function debug_exome_analysis()
Assign24MammalsGlobalConstants; AssignGeneralConstants; AssignStatsConstants; AssignRVASConstants;
X = load('temp_surface2.mat');

X.D.SFS.x_vec = X.x_vec_cell; X.D.SFS.p_vec = X.p_vec_cell; X.D.s_grid = X.s_vec; % overwrite SFS with UN-smoothed simulated results 
plot_params.figure_type = 1; plot_params.figs_dir = exome_data_figs_dir; 
plot_params.cum=1; plot_params.weighted = 1; plot_params.normalize=1;  plot_params.xlim = [10^(-6) 1]; % plot cumulative weighted allele frequency distribution
plot_allele_freq(X.s_vec, {X.D}, plot_params); 
%plot_allele_freq(X.s_vec([1 4:end]), {X.D}, plot_params); 

s_vec = [0 -logspace(-5, -1, 9)]; % set s-values for plotting and fitting
X.E = X.D; % create equilibrium model
X.E.init_pop_size_vec = {1000};
X.E.generations_vec = {2000};
X.E.expan_rate_vec = {1};
X.E.num_params_vec = [1 1 1]; X.E.num_stages=1;    X.E.index = 1;
X.E.iters = 250; 
X.E.name = 'equilibrium';
compute_flag = []; compute_flag.method = 'simulation'; compute_flag.smooth = 1;
[X.E.SFS.x_vec, X.E.SFS.p_vec, ...
    X.E.SFS.L, X.E.SFS.compute_time] = ...
    compute_allele_freq_spectrum_from_demographic_model( ...
    X.E, s_vec, compute_flag);
plot_params.xlim = [10^(-4) 1]; plot_allele_freq(s_vec, {X.E}, plot_params); title('SMOOTHED EQUILIBRIUM'); 
Y = load('temp_surface3.mat'); X.E.SFS.x_vec = Y.x_vec_cell; X.E.SFS.p_vec = Y.p_vec_cell; X.E.s_grid = Y.s_vec; % overwrite SFS with UN-smoothed simulated results 
plot_params.xlim = [10^(-4) 1]; plot_allele_freq(s_vec, {X.E}, plot_params);  title('UN-SMOOTHED EQUILIBRIUM');


exome_struct = get_exome_data_info(exome_data); % get metadata: file names, directories etc.
demography_file = [remove_suffix_from_file_name(exons_file) ...
    '_' 'AllPop' '_Demography.mat'];
demography_file = fullfile(spectrum_data_dir, ...
    exome_struct.data_str, 'AllPop', demography_file);

load(demography_file);
if(~isfield(Demographic_model{i_pop}, 'SFS')) % add SFS to demographic mode
    %    s_vec = [0 -logspace(-6, -2, 4)]; % light run - just for debugging
    Demographic_model{i_pop}.iters = 1000; % number of alleles to simulate !!
    Demographic_model{i_pop}.s_grid = [0 -logspace(-6, -2, 101)]; % s vector for interpolation
    compute_flag = []; compute_flag.method = 'simulation'; compute_flag.smooth = 1;
    [Demographic_model{i_pop}.SFS.x_vec, Demographic_model{i_pop}.SFS.p_vec, ...
        Demographic_model{i_pop}.SFS.L, SFS_compute_time] = ...
        compute_allele_freq_spectrum_from_demographic_model( ...
        Demographic_model{i_pop}, s_vec, compute_flag);
    save(demography_file, 'Demographic_model', 'max_LL_demographic_model');
end
