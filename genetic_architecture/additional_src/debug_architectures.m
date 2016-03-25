% Debug problem in architectures
% function X = debug_architectures()

% save('tmp_debug_arch.mat', 'architecture_str', ...
%         'f_vec', 'params_struct', 'p_x_vec', 'p_x_times_z_vec', 'iters', ...
%         'compute_method_flag', 'compute_marginals_flag'); 


% Compute mu for one architecture: 
p_vec = binopdf(0:6, 6, 0.2)

y = 0.5 *(1+tanh(8.23 .* ([0:6] - 0.281)))

y = tanh(8.23 .* ([0:6] - 0.281))

mumu2 = sum(p_vec .* y)

mumu = sum(p_vec .* y)
%mu^12


% % % load('tmp_debug_arch.mat');
% % % i=4;
% % % % [x_vec p_x_vec x_ind_vec x_ind_mat] = ...
% % % %     initilize_x_vec_constants(good_architectures(i).N, 0, good_architectures(i).f_vec(1));
% % % compute_marginals_flag = 1; 
% % % p_x_vec = [];  z = 1;
% % % p_x_vec = repmat(p_x_vec, iters, 1); % duplicate to match iters
% % % p_x_times_z_vec = p_x_vec .* z; % temporary vector to save time
% % % good_architectures(i).params_struct.min_freq = 0.08;
% % % good_architectures(i).params_struct.max_freq = 0.66;
% % % [v_marginal lods_ratio_marginal p_z_x_marginal] = ...
% % %     compute_architecture_statistics_marginal('and-of-sigmoids', ...
% % %     good_architectures(i).f_vec,  good_architectures(i).params_struct, p_x_vec, p_x_times_z_vec, 111, ...
% % %     'analytic', compute_marginals_flag); % compute several moments and other stuff for architecture






% Test the MZ-DZ twins to heritability estimate: 
XXX = load('temp_good_architecture.mat'); 

p_x_vec = []; p_x_times_z_vec = []; 
params_struct = XXX.good_architectures(2).params_struct
params_struct.min_freq = params_struct.min_freq(29);
params_struct.max_freq = params_struct.max_freq(29);

[mu V v_environment v_genetic mz_twin_risk H ] = ...
    compute_architecture_statistics('and-of-sigmoids', ...
    XXX.good_architectures(2).f_vec, params_struct, ...
    p_x_vec, p_x_times_z_vec, 100, ...
    'analytic')
[v_marginal GRR_marginal p_z_x_marginal h_add h_liability] = ...
    compute_architecture_statistics_marginal('and-of-sigmoids', ...
    XXX.good_architectures(2).f_vec, params_struct, p_x_vec, p_x_times_z_vec, 100, ...
    'analytic', mu)

[family_risk family_tree relative_risk sibs_genotype sibs_freqs ...
    H_from_twins h_add_from_twins h_liability_from_twins h_liability_from_twins_ADE] = ...
    compute_architecture_family_risk('and-of-sigmoids', ...
    XXX.good_architectures(2).f_vec, params_struct, 1000, ...
    'sampling', 2);
lambda_s = family_risk(end-1) / mu  % compute lambda_s
lambda_mz =  family_risk(end)/ mu

% twin_concordance_to_heritability % Compute heritability: 
h_liability_again = twin_concordance_to_heritability(lambda_mz, lambda_s, mu, 'ACE')
h_liability_based_on_MZ = familial_risk_to_heritability(lambda_mz, 'liability', mu, 1)
h_liability_based_on_DZ = familial_risk_to_heritability(lambda_s, 'liability', mu, 0.5)
lambda_dz_matching = heritability_to_familial_risk(h_liability_based_on_MZ, 'liability', mu, 0.5)
h_liability_should_match = twin_concordance_to_heritability(lambda_mz, lambda_dz_matching, mu, 'ACE')





h_liability_again_ADE = twin_concordance_to_heritability(lambda_mz, lambda_s, mu, 'ADE')

