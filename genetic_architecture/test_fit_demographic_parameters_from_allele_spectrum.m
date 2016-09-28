% Simulate spectrum from data. Then run to estimate demography. Check if
% correct demography recieves a good score (we could have non-unique solution) 
% function [D, max_LL] = test_fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec)

% Set demography
% New! create Demography structure
N = 500; 
num_generations = 100; 
expansion_factor = 1.02;
D.init_pop_size = [N -1]; %   N*2];
D.generations = [num_generations 1*num_generations]; %  num_generations]; % [num_generations 10*num_generations];
D.expan_rate = [expansion_factor 1]; %  0.99]; % [expansion_factor 1.0];
D.index = 1; % only one demography
D.add_new_alleles = 0; % use old simulation (track only alleles at start) - to be changed !!
n_sample = 200; % take a smaller sample size (can be equal to final population size)

N_vec = demographic_parameters_to_n_vec(D, 1);



% Simulate data: 
[x_vec, f_vec, k_vec, n_vec]  = compute_allele_freq_spectrum_from_demographic_model(D, 0, 'simulation', n_sample); % simulate from neutral model 


% now run and see if we get correct demography back 


[D_hat, max_LL] = fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec)


N_vec_hat = demographic_parameters_to_n_vec(D_hat, 1); 

figure; hold on; 
plot(N_vec); plot(N_vec_hat, 'r'); 
xlabel('Time (generations)'); ylabel('Population size'); 
legend({'True', 'Fitted'}); 
