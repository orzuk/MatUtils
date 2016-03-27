% Fit a simple demographic model to population allele frequency spectrum 
% We assume that observed alleles are not under selection (neutral) 
% 
% Input: 
% k_vec - number of dervied allele carriers observed in population for each allele
% n_vec - number of individuals profiled in population for each allele 
% 
% Output: 
% D - demographic model (expansion, population size etc.) 
% max_LL - loglikelihood of data for spectrum 
% 
function [D, max_LL] = fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec)

% We fit expansion model from allele-frequency data, assuming that all
% alleles are neutral
s = 0; % Assume no selection
alpha=0; % no mixture. Everything is neutral 
beta=[];

X = [k_vec n_vec]; % Represent counts in a packed form 

% Set different demographic models: 
D.bottleneck_size = logspace(1, 5, 50); % size of bottleneck
D.bottleneck_generations = logspace(1, 5, 50); % number of generations since bottleneck
D.init_pop_size = 2*10000; % ancestral population 
D.expansion_rate = logspace(0.001, 0.1, 50); % expansion rate per generation


for i=1:length(D) % Loop on D
    demographic_f_vec = compute_allele_freq_spectrum_from_demographic_model(D{i}); % Try a grid of different values
    
    % Compute likelihood. This is trivial one-class likelihood (no mixture bullshit) so should be fast !!!
    [log_like_mat(i)] = ... % compute likelihood (here vary only alpha)
        compute_two_class_log_likelihood(s, alpha, beta, rare_cumulative_per_gene, target_size_by_class_vec, N, ...
        X, [], [], [], [cur_is_null_vec is_null_mat(1:L_vec(i),i)], ...
        0, full_flag, num_individuals, D{i}); % don't include phenotype !!

end

[max_LL, max_J] = max(log_like_mat); % maximize likelihood 
D = D{max_J}; 



% Compute allele frequency distribution using Fisher-Wright model with changing population size. 
% Not sure if we need this here !!! 
% Input: 
% D - structure with demographic models
% s - selection coefficient 
% 
function [x_vec, f_vec] = compute_allele_freq_spectrum_from_demographic_model(D, s)

iters = 10^5; 
mu = 2*10^(-8); % set mutation rate 
init_str = 'newly_born'; 
compute_mode = 'simulation'
num_bins = 100; 


N_vec = demographic_parameters_to_n_vec(D); % compute population size at each generation

[freq_struct, absorption_struct, simulation_struct, N_vec, simulation_time] = ... % New: separate output to different structures
    FisherWrightSimulation(N_vec, mu, s, num_generations, expansion_factor, init_str, iters, compute_mode, num_bins); 

x_vec = freq_struct.x_vec; 
f_vec = p_vec; 


