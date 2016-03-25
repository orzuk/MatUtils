% Compute the phantom heritability for different generalized LP models
% function run_generalized_LP_different_P_norm_models()
figs_dir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/';

h_pathway = 0.8; % individual heritabilities 
N = 3; % number of pathways 
p_vec = [0.00001:0.05:9.911]; %  20 30 40 50]; 
iters = 1000000; % must be an integer produc of block size 
isoheritability_flag=0;

for i=1:length(p_vec)
    
[~, ~, ~, ~, ~, ...
    h_all_vec(i), h_pop_vec(i)] = ... %  ~, qtl_R_vec{i}] = ...
    compute_k_of_N_gaussian_statistics([0 0 0], [1 1 1], h_pathway, 0, iters, ...
    'EXP-SUM', p_vec(i), ... 'LP', p_vec(i), ...
    N, [], 'simulation', isoheritability_flag, {'ACE'});
end

[~, ~, ~, ~, ~, ...
    h_all_LP, h_pop_LP] = ... %  ~, qtl_R_vec{i}] = ...
    compute_k_of_N_gaussian_statistics([0 0 0], [1 1 1], h_pathway, 0, iters, ...
    'MAX', [], ...
    N, [], 'numeric', isoheritability_flag, {'ACE'}); % compute at L_infinity
[~, ~, ~, ~, ~, ...
    h_all_LP_sim, h_pop_LP_sim] = ... %  ~, qtl_R_vec{i}] = ...
    compute_k_of_N_gaussian_statistics([0 0 0], [1 1 1], h_pathway, 0, iters, ...
    'MAX', [], ...
    N, [], 'simulation', isoheritability_flag, {'ACE'}); % compute at L_infinity

h_pop_vec = cell2vec(h_pop_vec);
pi_phantom_vec = 1 - h_all_vec ./ h_pop_vec;
pi_phantom_LP = 1 - h_all_LP / h_pop_LP{1};

figure; hold on; w=5;
h_pop_vec_smoothed = smooth(p_vec, h_pop_vec, w)
h_all_vec_smoothed = smooth(p_vec, h_all_vec, w)
pi_phantom_vec_smoothed = smooth(p_vec, pi_phantom_vec, w)

plot_type = 'smoothed';
switch plot_type
    case 'smoothed'
        plot(p_vec, h_pop_vec_smoothed, 'linewidth', 2);
        plot(p_vec, h_all_vec_smoothed, 'r', 'linewidth', 2);
        plot(p_vec, pi_phantom_vec_smoothed, 'g', 'linewidth', 2);
    otherwise
        plot(p_vec, h_pop_vec);
        plot(p_vec, h_all_vec, 'r');
        plot(p_vec, pi_phantom_vec, 'g');
end

plot(10, h_pop_LP{1}, '*', 'markersize', 6);
plot(10, h_all_LP, '*r', 'markersize', 6);
plot(10, pi_phantom_LP, '*g', 'markersize', 6);
% plot(p_vec(end), h_pop_LP_sim{1}, '+');
% plot(p_vec(end), h_all_LP_sim, '+r');


xlabel('\alpha'); ylabel('heritability');
legend({'h_{pop}^2', 'h_{all}^2', '\pi_{phantom}'},4);
legend('boxoff');
title([repmat(' ', 1, 80) 'supp-fig6'], 'fontsize', 16, 'fontweight', 'bold');

my_saveas(gcf, fullfile(figs_dir, 'generalized_LP_heritability'), {'fig', 'jpg', 'epsc'});






