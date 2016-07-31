% Fit a simple demographic model to population allele frequency spectrum 
% We assume that observed alleles are not under selection (neutral) 
% 
% Input: 
% k_vec - number of dervied allele carriers observed in population for each allele
% n_vec - number of individuals profiled in population for each allele 
% 
% Output: 
% D - demographic model (expansion, population size etc.) 
% max_LL - log-likelihood of data for spectrum 
% 
function [D, max_LL] = fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec)

% We fit expansion model from allele-frequency data, assuming that all alleles are neutral
s = 0; % Assume no selection (synonymous) 
alpha=0; % No mixture. Everything is neutral (synonymous)
beta=[];  % effect size on phenotype IRRELEVANT! (we use only genpotypes) 

X = [k_vec n_vec]; % Represent counts in a packed form 

% Set different demographic models: 
D.init_pop_size_vec{1} = 2*10000; % 1. ancestral population 
D.init_pop_size_vec{2} = round(logspace(1, 4, 20)); % 2. bottleneck size 
D.init_pop_size_vec{3} = -1; % 3. population after bottleneck
%D.num_params_vec(1) = length(D.init_pop_size_vec{2}); 

D.generations_vec{1} = 500; % 4. burnout time 
D.generations_vec{2} = 10;  % 5. number of generations since bottleneck
D.generations_vec{3} = round(logspace(1, 4, 25));  % 6. number of generations for expansion 
%D.num_params_vec(2) = length(D.generations_vec{3}); 

D.expan_rate_vec{1} = 1; % fixed ancestral population size 
D.expan_rate_vec{2} = 1; % fixed buttleneck size 
D.expan_rate_vec{3} = logspace(0.001, 0.1, 30); % 5. expansion rate per generation
%D.num_params_vec(3) = length(D.expan_rate_vec{3}); 

D.num_params_vec = [length_cell(D.init_pop_size_vec) length_cell(D.generations_vec) length_cell(D.expan_rate_vec)];
D.num_stages = length(D.init_pop_size_vec);

%D.num_bottleneck_size = length(D.bottleneck_size);
%D.num_bottleneck_size = length(D.bottleneck_size);
%D.num_init_pop_size = length(D.bottleneck_size);
%D.num_expansion_rate = length(D.expansion_rate);


D.num_params = prod(D.num_params_vec); % how many parameters to enumerate


log_like_mat = zeros(size(D)); 
for i=1:D.num_params 
    
    length(D) % Loop on D: enumerate on many different demographic parameters 
    
    D.index = i;
    N_vec = demographic_parameters_to_n_vec(D, D.index); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation
    % Filter first unreasonable models !!! 
    if(max(N_vec) > 10^(10)) % 10 billion
        continue;
    end
    if(length(N_vec) > 10^4) % 10000 generations
        continue;
    end
    
    [demographic_x_vec, demographic_f_vec] = compute_allele_freq_spectrum_from_demographic_model(D, 0); % Try a grid of different values
    
    % Compute likelihood. This is trivial one-class likelihood (no mixture bullshit) so should be fast !!!
    rare_cumulative_per_gene = []; % set dummy variables 
    target_size_by_class_vec = []; % ??? 
    [log_like_mat(i)] = ... % compute likelihood (here vary only alpha)
        compute_two_class_log_likelihood(s, alpha, beta, rare_cumulative_per_gene, target_size_by_class_vec, N, ...
        X, [], [], [], [cur_is_null_vec is_null_mat(1:L_vec(i),i)], ...
        0, full_flag, num_individuals, D{i}); % don't include phenotype !!
end

[max_LL, max_J] = max(log_like_mat); % maximize likelihood 
D = D{max_J}; 



% Compute allele frequency distribution using Fisher-Wright model with changing population size. 
% (Not sure if we need this here !!! )
% Input: 
% D - structure with demographic models
% i - index of the demographic model to use out of the structure 
% s - selection coefficient 
% 
% Output: 
% x_vec - vector of x values (allele frequencies) at each generation 
% f_vec - vector of their frequencies at each generation 
% 
function [x_vec, f_vec] = compute_allele_freq_spectrum_from_demographic_model(D, s)

iters = 10^3; % number of alleles to simulate (start low to save time. As we refine demography fitting we increase this number) 
mu = 2*10^(-8); % set mutation rate 
init_str = 'newly_born'; 
compute_mode = 'simulation'; % for general demography
num_bins = 100; % used for binning in Fisher Right simulation

% i_vec = myind2sub(D.num_params_vec, length(D.num_params_vec), i); % create vector of indices
% D.generations = zeros(1, D.num_stages); D.expan_rate = zeros(1, D.num_stages); D.init_pop_size = zeros(1, D.num_stages); 
% for j=1:D.num_stages
%     D.init_pop_size(j) = D.init_pop_size_vec{j}(i_vec(j));     
%     D.generations(j) = D.generations_vec{j}(i_vec(j+D.num_stages)); 
%     D.expan_rate(j) = D.expan_rate_vec{j}(i_vec(j+2*D.num_stages)); 
% end
N_vec = demographic_parameters_to_n_vec(D, D.index); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation

num_final_generations = length(N_vec)-1; % simulation at the end 
[freq_struct, absorption_struct, simulation_struct, N_vec, simulation_time] = ... % New: separate output to different structures
    FisherWrightSimulation([], D, mu, s, init_str, iters, compute_mode, num_bins); 

x_vec = freq_struct.x_vec; 
f_vec = freq_struct.p_vec; 


