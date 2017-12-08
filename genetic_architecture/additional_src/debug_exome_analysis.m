%function debug_exome_analysis()
Assign24MammalsGlobalConstants; AssignGeneralConstants; AssignStatsConstants; AssignRVASConstants;
X = load('temp_surface.Fitted.African.mat'); %X = load('temp_surface2.mat');
X.N_vec = demographic_parameters_to_n_vec(X.D, X.D.index); 


X.D.SFS.x_vec = X.x_vec_cell; X.D.SFS.p_vec = X.p_vec_cell; X.D.s_grid = X.s_vec; % overwrite SFS with UN-smoothed simulated results 
plot_params.figure_type = 1; plot_params.figs_dir = exome_data_figs_dir; 
plot_params.cum=1; plot_params.weighted = 1; plot_params.normalize=1;  plot_params.xlim = [10^(-6) 1]; % plot cumulative weighted allele frequency distribution
Y = load('temp_surface3.mat'); Y.smooth_params.knots = 10; Y.smooth_params.y_fit = [0 logspace(-6, -1, 101)]; 
num_s = length(X.s_vec); x_vec = unique( [X.D.SFS.x_vec{:}] ); num_x = length(x_vec); p_mat = zeros(num_s, num_x);
for i_s = 1:length(Y.s_vec)
    [~, I, J] = intersect(x_vec, Y.x_vec_cell{i_s});
    p_mat(i_s,I) = Y.p_vec_cell{i_s}(J);
end
Y.smooth_params.plot=0;
% new: fit with cell array 
X.D.SFS.y_vec = rude(length_cell(X.D.SFS.x_vec), X.D.s_grid);


SM.name = 'Fitted.African-smoothed'; SM.s_grid = Y.smooth_params.y_fit; % s_vec; 
plot_params.weighted = 1; plot_params.cum = 1;  plot_params.log=[1 0]; plot_params.xlim = [10^(-6) 1];
plot_allele_freq(X.s_vec, {X.D}, plot_params); % plot once without producing histogram 


num_bins = min(X.N_vec(end), 1000); % Make histogram of p_vec 
x_bins = [0 unique(round(logspace(0, log10(2*X.N_vec(end)), num_bins)))];
x_bins2 = 0.5 * (x_bins + [x_bins(2:end)-1 x_bins(end)]);
for i=1:num_s
    X.D.SFS.save_p_vec{i} = X.D.SFS.p_vec{i};
    X.D.SFS.save_x_vec{i} = X.D.SFS.x_vec{i};
    X.D.SFS.p_vec{i} = weighted_hist(double(X.D.SFS.x_vec{i}), X.D.SFS.p_vec{i}, x_bins);
    X.D.SFS.p_vec{i} = X.D.SFS.p_vec{i} ./ [1 diff(x_bins)];
    X.D.SFS.x_vec{i} = x_bins2;
end
s_mesh = abs(mat2vec(repmat(X.s_vec, length(x_bins), 1)')); % figure; semilogx(x_bins, x_vals); 
%[SM.SFS.x_vec, SM.s_grid, SM.SFS.p_vec] = fit_monotonic_surface([X.D.SFS.x_vec{:}], abs(X.D.SFS.y_vec), [X.D.SFS.p_vec{:}], Y.smooth_params);  % constraints
[SM.SFS.x_vec, SM.s_grid, SM.SFS.p_vec] = fit_monotonic_surface([X.D.SFS.x_vec{:}], s_mesh, [X.D.SFS.p_vec{:}], Y.smooth_params);  % constraints

%[SM.SFS.x_vec, SM.s_grid, SM.SFS.p_vec] = fit_monotonic_surface(x_vec, abs(double(Y.s_vec)), p_mat, Y.smooth_params);  % constraints
SM.name = 'Fitted.African-smoothed'; SM.s_grid = Y.smooth_params.y_fit; % s_vec; 
plot_params.weighted = 1; plot_params.cum = 1;  plot_params.log=[1 0]; plot_params.xlim = [10^(-6) 1]; plot_params.hist = 1; plot_params.new_fig=1; 
plot_allele_freq(X.s_vec, {X.D}, plot_params); plot_params.new_fig=0;
plot_allele_freq(X.s_vec, {SM}, plot_params); plot_params.new_fig=1; %title('SMOOTHED EQUILIBRIUM AGAIN');
plot_params.weighted = 0; plot_params.cum = 0; plot_params.log=[1 1]; % Plot density (unweighted)
plot_allele_freq(X.s_vec, {X.D}, plot_params); plot_params.new_fig=0;
plot_allele_freq(X.s_vec, {SM}, plot_params); plot_params.new_fig=1;  %title('SMOOTHED EQUILIBRIUM AGAIN');
plot_params.weighted = 1; plot_params.cum = 0; plot_params.log=[1 1]; % Plot density (unweighted)
plot_allele_freq(X.s_vec, {X.D}, plot_params);  plot_params.new_fig=0;
plot_allele_freq(X.s_vec, {SM}, plot_params); % title('SMOOTHED EQUILIBRIUM AGAIN');


figure; loglog(X.D.SFS.save_x_vec{1}, (X.D.SFS.save_p_vec{1})); hold on;
loglog(X.D.SFS.x_vec{1}, (X.D.SFS.p_vec{1}), 'r'); legend('simulated', 'histogram');



return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss = -[0 logspace(-6, -1, 11)]; dd = []; 
N=10000; dd.SFS.x_vec = (1:(2*N)) ./ (2*N); 

for i=1:length(ss)
%    dd.SFS.z_vec(i) =  phi_s_integral(1-1/(2*N), ss(i)*4*N, 0) - phi_s_integral(1/(2*N), ss(i)*4*N, 0);
    dd.SFS.z_vec(i) = absorption_time_by_selection(ss(i), 1, N, 1/(2*N), 1-1/(2*N));
    dd.SFS.p_vec{i} = allele_freq_spectrum(dd.SFS.x_vec, ss(i), N, 0, 'linear') ./ dd.SFS.z_vec(i);    
end
dd.s_grid = ss; dd.name = 'Debug'; 
%figure; semilogx(xx, gg); 
parpar = []; parpar.cum=0; parpar.log=[0 1]; parpar.figure_type = 1; %parpar.xlim = [10^(-6) 1]; 
plot_allele_freq(ss, {dd}, parpar); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot_allele_freq(X.s_vec([1 4:end]), {X.D}, plot_params); 
X = load('temp_surface2.mat');
[X.D.SFS.x_vec, X.D.SFS.p_vec, ...
    X.D.SFS.L, X.D.SFS.compute_time] = ...
    compute_allele_freq_spectrum_from_demographic_model( ...
    X.D, Y.s_vec, compute_flag); % run again fitted expansion model (takes ~15 minutes)
plot_allele_freq(X.s_vec, {X.D}, plot_params); 

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
Y = load('temp_surface3.mat'); 
plot_params.xlim = [10^(-4) 1]; plot_allele_freq(Y.s_vec, {X.E}, plot_params); title('SMOOTHED EQUILIBRIUM'); 
X.E.SFS.x_vec = Y.x_vec_cell; X.E.SFS.p_vec = Y.p_vec_cell; X.E.s_grid = Y.s_vec; % overwrite SFS with UN-smoothed simulated results 
plot_params.weighted = 1; plot_params.cum = 1;  plot_params.xlim = [10^(-4) 1]; plot_allele_freq(Y.s_vec, {X.E}, plot_params);  title('UN-SMOOTHED EQUILIBRIUM');

% Perform smoothing again (outside):
num_s = length(Y.s_vec); x_vec = unique( [Y.x_vec_cell{:}] ); num_x = length(x_vec); p_mat = zeros(num_s, num_x);
for i_s = 1:length(Y.s_vec)
    [~, I, J] = intersect(x_vec, Y.x_vec_cell{i_s});
    p_mat(i_s,I) = Y.p_vec_cell{i_s}(J);
end
 Z = load('temp_surface2.mat'); Y.smooth_params.y_fit = abs(Z.D.s_grid); 
% Y.smooth_params.y_fit = abs(X.E.s_grid); 
Y.smooth_params.y_fit = abs([0 -logspace(-6, -1, 101)]);
[SM.SFS.x_vec, SM.s_grid, SM.SFS.p_vec] = fit_monotonic_surface(x_vec, abs(Y.s_vec), p_mat, Y.smooth_params);  % constraints
SM.name = 'equilibrium-smoothed'; SM.s_grid = Y.smooth_params.y_fit; % s_vec; 
plot_params.weighted = 1; plot_params.cum = 1; plot_params.xlim = [10^(-4) 1];  plot_allele_freq(Y.s_vec, {SM}, plot_params);  title('SMOOTHED EQUILIBRIUM AGAIN');
%plot_params.weighted = 1; plot_params.cum = 0;  plot_allele_freq(s_vec, {SM}, plot_params);  title('SMOOTHED EQUILIBRIUM DENSITY AGAIN');

save_p_vec = SM.SFS.p_vec(10,:);

figure; plot(save_p_vec -  SM.SFS.p_vec(10,:), '*')
figure; plot(save_p_vec, SM.SFS.p_vec(10,:), '*')
figure; semilogx(cumsum(save_p_vec), 'b'); hold on; semilogx(cumsum(SM.SFS.p_vec(10,:)), 'r--'); 
figure; plot(SM.SFS.x_vec .* save_p_vec, 'b'); hold on; plot(SM.SFS.x_vec .* SM.SFS.p_vec(10,:), 'r--'); 
figure; plot(cumsum(SM.SFS.x_vec .* save_p_vec), 'b'); hold on; plot(cumsum(SM.SFS.x_vec .* SM.SFS.p_vec(10,:)), 'r--'); 


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





%%%%%%%%%%% Power calculations %%%%%%%%%%%
% % a_vec = [0.05 2.5*10^(-6)];
% % b_vec = [0.01 0.05 0.1 0.5 0.9 0.95 0.99]; 
% % nu_mat = zeros(length(a_vec), length(b_vec)); 
% % for iii = 1:2
% %     for jjj = 1:length(b_vec)
% %         nu_mat(iii, jjj) = two_types_errors_to_non_centrality_parameter(a_vec(iii), b_vec(jjj));
% %     end
% % end
% % R_nu_mat = [ ['b \ a' num2cell(a_vec)]'  [num2cell(b_vec') num2cell(nu_mat')]' ]'
% % R_nu_mat = num2str_cell(R_nu_mat, 5);
% % savecellfile(R_nu_mat, 'nu_mat_tab.txt')
% % sprintf(R_nu_mat)
