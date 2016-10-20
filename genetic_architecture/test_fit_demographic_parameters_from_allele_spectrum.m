% Simulate spectrum from data. Then run to estimate demography. Check if
% correct demography recieves a good score (we could have non-unique solution)
% function [D, max_LL] = test_fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec)

% Set demography
% New! create Demography structure
test_population_to_sample=0;
N = 500;
num_generations = 100;
expansion_factor = 1.02;
D.init_pop_size = [N -1]; %   N*2];
D.generations = [num_generations 1*num_generations]; %  num_generations]; % [num_generations 10*num_generations];
D.expan_rate = [expansion_factor 1]; %  0.99]; % [expansion_factor 1.0];
D.index = 1; % only one demography
D.add_new_alleles = 1; % use NEW! simulation (track all alleles: alleles at start and newly born alleles) - to be changed !!
n_sample = 200; % take a smaller sample size (can be equal to final population size)

N_vec = demographic_parameters_to_n_vec(D, 1);

mu = 2*10^(-8) * 100000; % take effective mutation rate in a region 
D.mu=mu; 

% Simulate data:
[x_vec, p_vec, L_correction_factor, ~, k_vec, n_vec, weights_vec]  = ...
    compute_allele_freq_spectrum_from_demographic_model(D, 0, 'simulation', n_sample, mu); % simulate from neutral model

% D_equi = D; D_equi.expan_rate(:) = 1; 
% [x_vec_equi, p_vec_equi, k_vec_equi, n_vec_equi]  = compute_allele_freq_spectrum_from_demographic_model(D_equi, 0, 'simulation', n_sample); % simulate from neutral model
% p_vec_equi_analytic = 1 ./ x_vec_equi; p_vec_equi_analytic(1) = 0; p_vec_equi_analytic(end) = 0; 
% p_vec_equi_analytic=p_vec_equi_analytic./sum(p_vec_equi_analytic); 
% sum(p_vec_equi_analytic .* (x_vec_equi./(2*N)) .* (1-x_vec_equi./(2*N)) )
% sum(p_vec_equi .* (x_vec_equi./(2*N)) .* (1-x_vec_equi./(2*N)) )
% sum(p_vec .* (x_vec./(2*N_vec(end))) .* (1-x_vec./(2*N_vec(end))) )

% now run and see if we get correct demography back


[D_hat, max_LL] = fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec, weights_vec, mu, L_correction_factor, D)


N_vec_hat = demographic_parameters_to_n_vec(D_hat, 1);

figure; hold on;
plot(N_vec); plot(N_vec_hat, 'r');
xlabel('Time (generations)'); ylabel('Population size');
legend({'True', 'Fitted'});






% Temp tesT:
if(test_population_to_sample)
    n_sample = 25;
    [x_vec, p_vec, ~, ~, k_vec, n_vec, weights_vec] = compute_allele_freq_spectrum_from_demographic_model(D, 0, [], n_sample);
    
    allele_freq_vec = hist_to_vals(x_vec, p_vec ./ min(p_vec(p_vec>0)));
    
    num_alleles = length(allele_freq_vec);
    k_vec2 = population_to_sample_allele_freq(allele_freq_vec, 2*N_vec(end-1), n_sample);
    n_vec2 = repmat(n_sample, num_alleles, 1);
    
    
    
    %k_vec = population_to_sample_allele_freq(f_vec, N, n_sample)
    [sample_x_vec, sample_p_vec] = population_to_sample_allele_freq_distribution(x_vec, p_vec, n_sample);
    [hhh, ccc] = hist(k_vec, 0:n_sample); simulated_p_vec = hhh./sum(hhh);
    [hhh2, ccc2] = hist(k_vec2, 0:n_sample); simulated_p_vec2 = hhh2./sum(hhh2); figure;
    figure; bar([simulated_p_vec', simulated_p_vec2', sample_p_vec']); legend('sampled', 'sampled2', 'average');
    
    sample_p_vec
    simulated_p_vec
    simulated_p_vec2
end
