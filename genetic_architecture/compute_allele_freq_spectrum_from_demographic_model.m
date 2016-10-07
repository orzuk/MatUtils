% Compute allele frequency distribution using Fisher-Wright model with changing population size. 
% (Not sure if we need this here !!! )
% Input: 
% D - structure with demographic models
% s - selection coefficient 
% compute_flag - 'simulation' (default) or 'moments' (computation based on moments)
% 
% Output: 
% x_vec - vector of x values (allele frequencies) at each generation 
% f_vec - vector of their frequencies at each generation 
% 
function [x_vec, f_vec, k_vec, n_vec] = compute_allele_freq_spectrum_from_demographic_model(D, s, compute_flag, n_sample)

if(~exist('compute_flag', 'var') || isempty(compute_flag))
    compute_flag = 'simulation'; 
%    compute_mode = 'simulation'; % for general demography
end

iters = 10^3; % number of alleles to simulate (start low to save time. As we refine demography fitting we increase this number) 
mu = 2*10^(-8); % set mutation rate 
init_str = 'newly_born'; 
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

switch compute_flag
    case {'simulation', 'simulations', 'numeric'}
        [freq_struct, absorption_struct, simulation_struct, N_vec, simulation_time] = ... % New: separate output to different structures
            FisherWrightSimulation([], D, mu, s, init_str, iters, compute_flag, num_bins);
        x_vec = freq_struct.x_vec{end-1}; % why don't take last one? 
        f_vec = freq_struct.p_vec{end-1};

        num_alleles = length(simulation_struct.q(:,end)); 
        k_vec = population_to_sample_allele_freq(simulation_struct.q(:,end), N_vec(end-1), n_sample);
        n_vec = repmat(n_sample, num_alleles, 1); 
        
    case 'moments'
        
        N_vec = demographic_parameters_to_n_vec(D, 1);
        [mu_vec_expansion_analytic] = FisherWright_Compute_SFS_Moments(N_vec, 0, max_k); % compute moments with Formulas from Ewens
    
        % Estimate density using the max-entropy method? NO! just compute moments !
        x_vec = (1:(2*N-1)) ./ (2*N); % use initial pop size for resolution 
        [lambda_max_ent] = ...  % Fit density using moments % , lambda0)
            me_dens2(mu_vec_expansion_analytic(2:end) ./ mu_vec_expansion_analytic(1), x_vec); % normalize by 0th moment
        
        f_vec = zeros(1, 2*N-1);
        for k=1:(max_k)
            f_vec = f_vec + lambda_max_ent(k) .* x_vec .^ (k-1);
        end
        f_vec = exp(-f_vec) ./ (x_vec .* (1-x_vec)); % Compute f. How to normalize?
        
end


