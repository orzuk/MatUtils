% Set vector of population sizes for different demographic models
%
% Input:
% demographic_str - string representing the demography
% num_generation - number of 'burn-in' generations to start with
%
% Output:
% N_vec - vector of population sizes
%
function N_vec = create_demographic_model(demographic_str, num_generations)


switch lower(demographic_str)
    case 'finland'
        
    case 'iceland'
        
    case 'europe'
        
    case 'equilibrium'
        N_vec = repmat(10000, num_generations, 1);
        
    case 'expansion1'
        
    case 'expansion2'
        
    case 'tennensen'
        
    case 'coventry'
        
    case 'nelson'
        
    case 'gazava' % main model
        N = 1000; % initial population size
        burn_in = 20; % 20 generations at equalibrium
        growth_gen = 400; % growth
        N_final = 10000000;
        
        N_vec = [repmat(N, 1, burn_in) start_end_to_n_vec_internal(N, N_final, growth_gen)]'; % set N 
end


function N_vec = start_end_to_n_vec_internal(N_init, N_final, growth_gen)

N_vec = logspace(log10(N_init), log10(N_final), growth_gen); 



