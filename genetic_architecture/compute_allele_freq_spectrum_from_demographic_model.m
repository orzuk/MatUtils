% Compute allele frequency distribution using Fisher-Wright model with changing population size.
% (Not sure if we need this here !!! )
% Input:
% D - structure with demographic models
% s - selection coefficient
% compute_flag - 'simulation' (default) or 'moments' (computation based on moments)
% n_sample - # of individuals in a sample 
% mu - regional mutation rate 
%
% Output:
% x_vec - vector of x values (allele frequencies) at each generation
% p_vec - vector of their frequencies at each generation
% k_vec - alleles in sample
% n_vec - sample sizes
%
function [x_vec, p_vec, k_vec, n_vec, L_correction_factor, compute_time] = ...
    compute_allele_freq_spectrum_from_demographic_model(D, s, compute_flag, n_sample, mu)

compute_time=cputime;
if(~exist('compute_flag', 'var') || isempty(compute_flag))
    compute_flag = 'simulation';
    %    compute_mode = 'simulation'; % for general demography
end
if(~exist('init_str', 'var') || isempty(init_str)) % default is start at equilibrium
    init_str = 'equilibrium';
end

if(~exist('mu', 'var') || isempty(mu))
    mu = 2*10^(-8); % set mutation rate
end

if(~isfield('iters', D))
%    D.iters = 5000;
    D.iters = 1000; % number of alleles to simulate (start low to save time. As we refine demography fitting we increase this number)
end
D.num_bins = 100; % used for binning in Fisher Right simulation
D.compute_absorb = 0; % no need for extra computation!!! 

N_vec = demographic_parameters_to_n_vec(D, D.index); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation
if(~exist('n_sample', 'var') || isempty(n_sample))
   n_sample =  2*N_vec(end-1);
end


num_final_generations = length(N_vec)-1; % simulation at the end
L_correction_factor=[]; 

switch compute_flag
    case {'simulation', 'simulations', 'numeric'}
        [freq_struct, ~, simulation_struct, N_vec, simulation_time] = ... % New: separate output to different structures
            FisherWrightSimulation([], D, mu, s, init_str, D.iters, compute_flag, D.num_bins);
        x_vec = freq_struct.x_vec{end-1}; % why don't take last one?
        p_vec = freq_struct.p_vec{end-1};
        L_correction_factor = simulation_struct.L_correction_factor;
                
        [sample_x_vec, sample_p_vec] = population_to_sample_allele_freq_distribution(x_vec, p_vec, n_sample); % compute distribution at sample (no further sampling here!)
        if(~isfield(simulation_struct, 'weights'))
            simulation_struct.weights = []; 
        end
        % currently: round p_vec to integers! (could be innaccurate, and also have many alleles for large N)
        allele_freq_vec = hist_to_vals(x_vec, round(p_vec .* simulation_struct.num_simulated_polymorphic_alleles_vec(end))); % compute alleles at sample ?
        num_alleles = length(allele_freq_vec);
        pop_to_sample_t = cputime; 
        k_vec = population_to_sample_allele_freq(allele_freq_vec, 2*N_vec(end-1), n_sample); % simulate a sample from population 
        pop_to_sample_t = cputime-pop_to_sample_t; 
        fprintf('Converted %d alleles to sample freq. time=%f\n', num_alleles, pop_to_sample_t); 
        
        n_vec = repmat(n_sample, num_alleles, 1);
        
    case 'moments'
        
        N_vec = demographic_parameters_to_n_vec(D, 1);
        [mu_vec_expansion_analytic] = FisherWright_Compute_SFS_Moments(N_vec, 0, max_k); % compute moments with Formulas from Ewens
        
        % Estimate density using the max-entropy method? NO! just compute moments !
        x_vec = (1:(2*N-1)) ./ (2*N); % use initial pop size for resolution
        [lambda_max_ent] = ...  % Fit density using moments % , lambda0)
            me_dens2(mu_vec_expansion_analytic(2:end) ./ mu_vec_expansion_analytic(1), x_vec); % normalize by 0th moment
        
        p_vec = zeros(1, 2*N-1);
        for k=1:(max_k)
            p_vec = p_vec + lambda_max_ent(k) .* x_vec .^ (k-1);
        end
        p_vec = exp(-p_vec) ./ (x_vec .* (1-x_vec)); % Compute f. How to normalize?
        
end

compute_time=cputime-compute_time;
