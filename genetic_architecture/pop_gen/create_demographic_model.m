% Set vector of population sizes for different demographic models
%
% Input:
% demographic_str - string representing the demography
% num_generation - number of 'burn-in' generations to start with
%
% Output:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              % Output:
% N_vec - vector of population sizes
% D - structure of demographic models 
% 
function [N_vec, D] = create_demographic_model(demographic_str, num_generations)

if(~exist('num_generations', 'var') || isempty(num_generations)) % set default 'burn-in' 
    num_generations = 10;
end
N = 10000; % default equilibrium pop. size 
D.index = 1; % only one demography
D.add_new_alleles = 1; % Default: add new alleles each generation (this is part of demography?) % use old simulation (track only alleles at start) - to be changed !!

switch lower(demographic_str)
    case 'two-stage-expan'
        N = 500; 
        D.init_pop_size = [N -1];
        D.generations = [num_generations num_generations]; 
        D.expan_rate = [1.02 1]; 
    
    case 'small-expan'
        N = 200; 
        D.init_pop_size = [N -1];
        D.generations = [num_generations 100]; 
        D.expan_rate = [1 1.01];
    case 'finland'
        D.init_pop_size = [N -1 50];
        D.generations = [num_generations 900 100]; 
        D.expan_rate = [1 1.00462 1.122];
    case 'iceland'
        D.init_pop_size = [N -1 500];
        D.generations = [num_generations 950 50]; 
        D.expan_rate = [1 1.00462 1.132];        
    case 'europe'
        D.init_pop_size = [N 775 -1];
        D.generations = [num_generations 1230 50]; 
        D.expan_rate = [1 1.00396 1.094];
    case 'equilibrium'
        D.init_pop_size = N;
        D.generations = num_generations;
        D.expan_rate = 1;        
    case 'small-equilibrium'
        N = 500; 
        D.init_pop_size = N;
        D.generations = num_generations;
        D.expan_rate = 1;    
    case 'tiny-equilibrium'
        N = 5; 
        D.init_pop_size = N;
        D.generations = num_generations;
        D.expan_rate = 1;    
    case 'expansion1'
        D.init_pop_size = [N -1 -1];
        D.generations = [num_generations 1000]; 
        D.expan_rate = [1 1.00462];
    case 'expansion2'
        D.init_pop_size = [N -1 -1];
        D.generations = [num_generations 1230 50]; 
        D.expan_rate = [1 1.00197 1.094];

    case 'tennensen'
        
    case 'coventry'
        
    case 'nelson'
        
    case 'gazava' % main model
        N = 1000; N_final = 10000000;
        D.init_pop_size = [N -1]; % initial population size
        D.generations = [20 400]; % 20 burn-in, 400 growth
        D.expan_rate = [1 (N_final/N)^(1/(400-1))];
%        N_vec = [repmat(N, 1, burn_in) start_end_to_n_vec_internal(N, N_final, growth_gen)]'; % set N 
end
D.num_stages = length(D.generations);
% Finaly convert model to N vec
N_vec = demographic_parameters_to_n_vec(D, 1); 
D.total_generations = length(N_vec)-1;
D.name = demographic_str; 

% Internal function: interpolate population size between initial and final sizes 
%function N_vec = start_end_to_n_vec_internal(N_init, N_final, growth_gen)

%N_vec = logspace(log10(N_init), log10(N_final), growth_gen); 



