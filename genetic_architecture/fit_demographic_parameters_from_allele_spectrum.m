% Fit a simple demographic model to population allele frequency spectrum 
% We assume that observed alleles are not under selection (neutral) 
% 
% Input: 
% k_vec - number of dervied allele carriers observed in population for each allele
% n_vec - number of individuals profiled in population for each allele 
% 
% Output: 
% D - demographic model (expansion, population size etc.) 
% max_LL - loglikelihood of data for 
% 
function [D max_LL] = fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec)

% We fit expansion model from allele-frequency data, assuming that all
% alleles are neutral
s = 0; % Assume no selection
alpha=0; % no mixture. Everything is neutral 
beta=[];

X = [k_vec n_vec]; % Represent counts in a packed form 

for i=1:length(D) % Loop on D
    demographic_f_vec = compute_allele_freq_spectrum_from_demographic_model(D{i}); % Try a grid of different values
    
    % Compute likelihood. This is trivial one-class likelihood (no mixture bullshit) so should be fast !!!
    [log_like_mat(i)] = ... % compute likelihood (here vary only alpha)
        compute_two_class_log_likelihood(s, alpha, beta, rare_cumulative_per_gene, target_size_by_class_vec, N, ...
        X, [], [], [], [cur_is_null_vec is_null_mat(1:L_vec(i),i)], ...
        0, full_flag, num_individuals, D{i}); % don't include phenotype !!

end

[max_LL max_J] = max(log_like_mat); % maximize likelihood 
D = D{max_J}; 



% Compute allele frequency using Fisher-Wright model with changing
% population size. 
% Not sure if we need this here !!! 
function compute_allele_freq_spectrum_from_demographic_model(D)

[q x_vec p_vec p_vec_analytic het_vec N_vec ...
    absorb_time_vec fixation_time_vec loss_time_vec total_het_vec ...
    frac_polymorphic_vec prob_fixation frac_alleles_kept_vec frac_het_kept_vec prob_site_polymorphic ...
    all_new_x_vec all_new_p_vec all_new_het_vec num_alleles_vec simulation_time] = ...
    FisherWrightSimulation(N, mu, s, num_generations, expansion_factor, init_str, iters, compute_mode, num_bins)


