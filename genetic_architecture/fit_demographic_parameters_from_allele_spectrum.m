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
D.init_pop_size_vec{2} = round(logspace(1, 5, 20)); % 2. bottleneck size 
D.init_pop_size_vec{3} = -1; % 3. population after bottleneck
%D.num_params_vec(1) = length(D.init_pop_size_vec{2}); 

D.generations_vec{1} = 500; % 4. burnout time 
D.generations_vec{2} = 10;  % 5. number of generations since bottleneck
D.generations_vec{3} = logspace(1, 5, 50);  % 6. number of generations for expansion 
%D.num_params_vec(2) = length(D.generations_vec{3}); 

D.expan_rate_vec{1} = 1; % fixed ancestral population size 
D.expan_rate_vec{2} = 1; % fixed buttleneck size 
D.expan_rate_vec{3} = logspace(0.001, 0.1, 50); % 5. expansion rate per generation
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
    [demographic_x_vec, demographic_f_vec] = compute_allele_freq_spectrum_from_demographic_model(D, i, 0); % Try a grid of different values
    
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
function [x_vec, f_vec] = compute_allele_freq_spectrum_from_demographic_model(D, i, s)

iters = 10^5; % number of alleles to simulate 
mu = 2*10^(-8); % set mutation rate 
init_str = 'newly_born'; 
compute_mode = 'simulation'; % for general demography
num_bins = 100; % used for binning in Fisher Right simulation

i_vec = myind2sub(D.num_params_vec, length(D.num_params_vec), i); % create vector of indices
D.generations = zeros(1, D.num_stages); D.expan_rate = zeros(1, D.num_stages); D.init_pop_size = zeros(1, D.num_stages); 
for j=1:D.num_stages
    D.generations(j) = D.generations_vec{j}(i_vec(j)); 
    D.expan_rate(j) = D.expan_rate_vec{j}(i_vec(j+D.num_stages)); 
    D.init_pop_size(j) = D.init_pop_size_vec{j}(i_vec(j+2*D.num_stages));     
end

N_vec = demographic_parameters_to_n_vec(D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation

num_final_generations = length(N_vec)-1; % simulation at the end 
[freq_struct, absorption_struct, simulation_struct, N_vec, simulation_time] = ... % New: separate output to different structures
    FisherWrightSimulation(N_vec, mu, s, num_final_generations, [], init_str, iters, compute_mode, num_bins); 

x_vec = freq_struct.x_vec; 
f_vec = freq_struct.p_vec; 


% Convert a set of demographic parameters into a vector of population sizes
% Input: 
% D - a structure with demographic parameters 
% i - index stating which demographic parameters to take  
%
% Output: 
% N_vec - a vector of population size at each generation 
% 
function N_vec = demographic_parameters_to_n_vec(generations_vec, expan_rate_vec, init_pop_size_vec) %   (D, i)

% Here we expand: gen_vec, expan_rate_vec, init_pop_size_vec  to one vector N_vec containing population size at each generation
num_generations = sum(generations_vec); 

N_vec = zeros(num_generations, 1); 

ctr=0;
for i=1:length(generations_vec) % loop on #generation
   if(init_pop_size_vec(i) == -1)
       init_size = N_vec(ctr); 
   else
       init_size = init_pop_size_vec(i); 
   end
   N_vec((ctr+1):(ctr+generations_vec(i))) = ceil(init_size * expan_rate_vec(i) .^ (1:generations_vec(i))); % get (rounded) population size 
   ctr = ctr + generations_vec(i); 
end

% % % i_vec = vec_to_array_ind(i, D.num_parmas); % create vector of indices

% % % num_generations = D.init_pop_generations(i_vec(2)) + D.bottleneck_generations(i_vec(4)) + D.expansion_generations(i_vec(6)); 
% % % N_vec = zeros(D.init_pop_generations(
% % % 
% % % num_generations(i_vec(1)), 1); % set number of generations 
% % % N_vec(1:D.init_pop_generations(i_vec(3))) = D.init_pop_size(i_vec(2)); % starting population size 
% % % 
% % % %i_gen=1; 
% % % 
% % % N_vec((D.init_pop_generations(i_vec(3))+1):D.init_pop_generations(i_vec(3))+D.bottleneck_generations(i_vec(5))) = D.bottleneck_size(i_vec(4)); 
% % % i_gen = D.init_pop_generations(i_vec(3))+D.bottleneck_generations(i_vec(5)) + 1; 
% % % 
% % % for j=1:length(D.expansion_generations(i_vec(7))) % loop on changes in population size rate     
% % %     i_gen = i_gen + D.bottleneck_generations(i_vec(4));
% % %     N_vec(i_gen) = round(N_vec(i_gen-1) * D.expansion_rate(i_vec(6))); i_gen = i_gen+1; % set population size     
% % % end

%%%D.bottleneck_size = logspace(1, 5, 50); % size of bottleneck
%%%D.bottleneck_generations = logspace(1, 5, 50); % number of generations since bottleneck
%%%D.init_pop_size = 2*10000; % ancestral population 
%%%D.expansion_rate = logspace(0.001, 0.1, 50); % expansion rate per generation


