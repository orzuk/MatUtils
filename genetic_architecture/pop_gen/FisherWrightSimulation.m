% Simulate a generalized Fisher-Wright model with selection and mutations.
% We allow expansion, and keep track of both old and new distribution of allele frequencies
%
% Input:
% region_flag - NEW! flag saying if to simulate a region with fixed size
%               (and then simulate iters such regions!) or simulate a fixed number of alleles (default)
% D - NEW! Demographic model. Contains population size at each generation !!!
% mu - mutation rate (per site/region per generation)
% s - selection coefficient. NOTE: Should be NEGATIVE for deleterious alleles
% init_str - initialization of allele frequencies (equilibrium or newly-born alleles)
% iters - how many times to perform simulation
% compute_mode - how to compute (simulation/markov-chain numeric computation/analytic)
% num_bins - # of bins in returend histograms (when N is large we want a coarse-grained version of allele frequencies) - NOT USED YET!!!
% plot_flag - if plotting SFS
%
% Output:
% freq_struct - structure with output on allele frequencies. Fields:
%               x_vec - vector of allele frequencies for each generation
%               p_vec - vector of their probabilities (NOT counts: how many
%               alleles present) for each generation . QU: should be counts or histogram?
%               p_vec_equilibrium_analytic - vector of allele frequencies computed analytically
%               het_vec - vector of heterozygosities
% absorption_struct - structure with output on absorption times. Fields:
%               absorption_time_given_init_freq_vec - distribution of
%               absorption times conditioned on initial frequencies
%               fixation_time_given_init_freq_vec - distribution of
%               fixation times conditioned on initial frequencies (and fixation)
%               loss_time_given_init_freq_vec - distribution of time to
%               loss conditioned on initial frequncies (and loss)
%               total_het_at_each_generation_vec - total heterozygosity present in population in each generation
% simulation_struct - structure with output on more simulation details. Fields:
%               frac_polymorphic_vec - fraction of alleles which are kept
%               polymorphic at each generation from starting alleles
%               prob_fixation - probability that an allele which reached absorption is fixed (compared to lost)
%               frac_old_alleles_survived_vec - fraction of old alleles which are still polymorphic at each generation
%               frac_het_kept_vec - fraction of old heterozygosity still present (at equilibrium)
%               prob_site_polymorphic_at_equilibrium - probability that a random site in the genome is polymorphic,
%               at the start of the model (equilibrium, before expansion)
%               all_new_x_vec - allele frequency values for newly born alleles (like x_vec, but just for new alleles)
%               all_new_p_vec - probability distribution for newly born alleles (like p_vec, but just for new alleles)
%               all_new_het_vec - heterozygosity distribution for newly born alleles
%               all_old_x_vec - allele frequency values for old alleles (like x_vec, but just for old alleles)
%               all_old_p_vec - probability distribution for old alleles (like p_vec, but just for old alleles)
%               all_old_het_vec - heterozygosity distribution for old alleles
%               num_simulated_polymorphic_alleles_vec - how many polymorphic alleles per generation were simulated (when we choose the simulated flag)
% N_vec - population size at each generation
% simulation_time - how much did the time computation took
%
function [freq_struct, absorption_struct, simulation_struct, N_vec, simulation_time] = ... % New: separate output to different structures
    FisherWrightSimulation(region_flag, D, mu, s, init_str, iters, compute_mode, num_bins, plot_flag)

absorption_time_given_init_freq_vec = []; q = [];
simulation_time =  cputime;
compute_matrix = 0; % compute transition matrix, absorption time etc.
if(~exist('region_flag', 'var') || isempty(region_flag))
    region_flag = 0; % % default: simulated a predefined number of alleles
end
if(~exist('init_str', 'var') || isempty(init_str))
    init_str = 'equilibrium'; % 'newly_born'; % default: start at equilibrium
end
two_side_flag = 0; % use derived allele frequency
if(~exist('plot_flag', 'var') || isempty(plot_flag))
    plot_flag = 0; % default: don't plot anything
end
if(~isfield(D, 'index')) % assume only one demography
    D.index = 1;
end
if(~isfield(D, 'add_new_alleles')) % generate new alleles in each generation
    D.add_new_alleles = 1;
end

% Extract all infromation from Demographic model
N_vec = demographic_parameters_to_n_vec(D, D.index);  % compute population size at each generation
N = N_vec(1); num_generations = length(N_vec)-1; max_N = max(N_vec);
prob_site_polymorphic_at_equilibrium = (2*N*mu) * 2 * absorption_time_by_selection(s, 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start
mu_vec_analytic = [];  mu_vec_equilibrium = [];

%heterozygosity_per_site = 4*N*mu * absorption_time_by_selection(s, 1, N, 0, 1, 'var');  % heterozygosity per site
%p_vec = cell(num_generations+1, 1); %total_het_at_each_generation_vec = zeros(num_generations,1);
num_moments = min(floor(N_vec(1)/3), 5);
switch compute_mode % how to compute 
    case {'analytic', 'moment', 'moments'} % compute solution based on Jacobi polynomials (See e.g. Kryukov et al. PNAS 2009).
        [freq_struct, absorption_struct, simulation_struct, ...
            num_absorptions ,num_fixations, ...
            num_absorptions_by_generation_vec, count_vec, num_losses] = ...
            compute_analytic_forward_FisherWright_internal( ...
            s, mu, two_side_flag, num_generations, N_vec, init_str, num_moments);

    case 'numeric' % here compute everything by matrix multiplications
        [freq_struct, absorption_struct, simulation_struct, ...
            num_absorptions, num_fixations, ...
            num_absorptions_by_generation_vec, count_vec, ...
            num_losses, M] = ...
            compute_numeric_forward_FisherWright_internal( ...
            s, mu, two_side_flag, N_vec, init_str, compute_matrix, num_bins);


    case 'simulation' % should be a sub-routine
        [freq_struct, absorption_struct, simulation_struct, q, ...
            num_absorptions, num_fixations, num_losses, ...
            num_absorptions_by_generation_vec, num_fixations_by_generation_vec, num_losses_by_generation_vec, ...
            count_vec] = simulate_forward_FisherWright_internal(D, s, mu, two_side_flag, iters, init_str);
        %        num_effective_iters = sum(p_vec{end}(2:end-1)) % total number of single-generation single-allele steps performed (?)
        
end % switch compute mode

% Fill additional statistics
het_vec = freq_struct.p_vec; % initilize distributions
moments_mat = zeros(num_generations, num_moments); het_moments_mat = moments_mat;
for j=1:num_generations % heterozygosity vector
    het_vec{j} = 2 .* freq_struct.p_vec{j} .* vec2row(freq_struct.x_vec{j}./(2*N_vec(j)) .* (1-freq_struct.x_vec{j}./(2*N_vec(j))));
    % NEW! compute moments!!!!
    central_flag = 0; % take non-central moments
    for k=1:num_moments
        moments_mat(j,k) = moment_hist(freq_struct.x_vec{j} ./ (2*N_vec(j)), freq_struct.p_vec{j}, k, central_flag); % compute moments of SFS distribution. Start from 1st moment!!
        het_moments_mat(j,k) = moment_hist(freq_struct.x_vec{j} ./ (2*N_vec(j)), het_vec{j}, k-1, central_flag, 0); % compute (k-1)-th (UN-Normalized!!!) moments of heterozygosity distribution. Start from zero-th moment!!
    end
end
final_x_vec = (1:2*N_vec(num_generations)-1) ./ (2*N_vec(num_generations)); % set new coordinates
freq_struct.p_vec_equilibrium_analytic = exp( allele_freq_spectrum([0 final_x_vec 1], s, N, two_side_flag, 'log') ); % compute analytic approxiamtion (valid only for constant population size)
freq_struct.p_vec_equilibrium_analytic = normalize_hist(final_x_vec, freq_struct.p_vec_equilibrium_analytic(2:end-1)); % normalized

all_new_p_vec = []; all_new_x_vec = [];
switch init_str % unite distributions into one
    case 'newly_born' % we start at newly born alleles
        for j=1:num_generations
            [all_new_x_vec, all_new_p_vec] = union_with_counts(all_new_x_vec, all_new_p_vec, ...
                freq_struct.x_vec{j}, freq_struct.p_vec{j} .* 2*mu * N_vec(num_generations-j+1)); % compute weighted sum
        end
        [all_new_x_vec, sort_perm] = sort(all_new_x_vec);
        all_new_p_vec = all_new_p_vec(sort_perm);
        all_new_het_vec = 2.* all_new_p_vec .* ...
            all_new_x_vec./(2*N_vec(num_generations)) .* (1-all_new_x_vec/(2*N_vec(num_generations)));
    otherwise % we start at equilibrium
        all_new_het_vec = [];
        all_old_x_vec = freq_struct.x_vec{num_generations} .* (2*N_vec(num_generations));
        all_old_p_vec = freq_struct.p_vec{num_generations};
        all_old_het_vec = het_vec{num_generations};
end
absorption_struct.prob_fixation = num_fixations / num_absorptions;
frac_old_alleles_survived_vec = [1 vec2row(cumprod(simulation_struct.frac_polymorphic_vec))]; % take total # of alleles left
frac_het_vec = 99999; % this is the fraction of heterozygosity that is retained at each generation (TEMP!!!)
absorption_struct.frac_het_kept_vec = [1 vec2row(cumprod(frac_het_vec))]; % take total # of alleles left
if(plot_flag) % Plot individual trajectories
    plot_expansion_internal(N, mu, s, q, freq_struct.x_vec, freq_struct.p_vec, freq_struct.p_vec_equilibrium_analytic);
end % if plot

% Build output structures
% New! only need to add a few fields !!! 
% % % % % % freq_struct = var2struct(x_vec, p_vec, p_vec_equilibrium_analytic, ...
% % % % % %     final_x_vec, het_vec, total_het_at_each_generation_vec, ...
% % % % % %     all_new_x_vec, all_new_p_vec, all_new_het_vec, compute_mode, ...
% % % % % %     prob_site_polymorphic_at_equilibrium, prob_site_polymorphic_at_end, ...
% % % % % %     mu_vec_analytic, mu_vec_equilibrium, moments_mat, het_moments_mat);
% % % % % % absorption_struct = var2struct(absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, loss_time_given_init_freq_vec, ...
% % % % % %     frac_polymorphic_vec, prob_fixation, frac_old_alleles_survived_vec, frac_het_kept_vec);
% % % % % % switch compute_mode % add simulations details
% % % % % %     case 'simulation'
% % % % % %         simulation_struct = var2struct(q, num_simulated_polymorphic_alleles_vec, weights, ...
% % % % % %             frac_polymorphic_vec, prob_fixation, frac_het_kept_vec, L_correction_factor); % add many things here
% % % % % %         
% % % % % %     otherwise % do not hold simulation information
% % % % % %         simulation_struct = [];
% % % % % %         simulation_struct.L_correction_factor = L_correction_factor;
% % % % % %         simulation_struct.num_simulated_polymorphic_alleles_vec = 10000; % large to minimize rounding error
% % % % % % end
simulation_time = cputime - simulation_time; % compute total time



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%% Internal function for simulation
% Internal function for simulation of alleles to get frequency distribution.
% Simulate until we've got enough 'good alleles'
%
% Input:
% N_vec - population size at each generation (number of INDIVIDUALS!)
% s - selection coefficient. Should be NEGATIVE for deleterious alleles
% mu - mutation rate: per REGION per individual per generation
% two_side_flag - derived/minor allele frequency
% iters - how many alleles to output
% num_generations - total number of generations to simulate
% init_str - how to start (equilibrium or new allele)
%
% Output:
% q - matrix with #of carriers simulated. Dimension: iters X num-generations
% weights - weight of different generations for simulation
% p_vec - vector of probabilities at each generation
% x_vec - vector of num. carriers at each generation
% total_het_at_each_generation_vec - total heterozygosity at each generation
% num_absorptions - number of absorptions (old alleles which dissappeared)
% num_fixations - number of fixations (old alleles which fixate in the population)
% num_absorptions_by_generation_vec - number of alleles which were absorped per-generation
% fixation_time_given_init_freq_vec - time to fixation for each allele which fixed
% num_losses - number of losses (old alleles which are lost in the population)
% loss_time_given_init_freq_vec - at what generation each loss occurs
% num_simulated_polymorphic_alleles_vec -
% absorption_time_given_init_freq_vec - time to absorptions for each allele which fixed
% count_vec - vector counting the .. ???
%
% [q, weights, x_vec, p_vec, total_het_at_each_generation_vec, ...
%     num_absorptions, num_fixations, num_losses, ...
%     num_absorptions_by_generation_vec, num_fixations_by_generation_vec, num_losses_by_generation_vec, ...
%     absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, loss_time_given_init_freq_vec, ...
%     num_simulated_polymorphic_alleles_vec, count_vec, prob_site_polymorphic_at_end, L_correction_factor] = ...
%     simulate_forward_FisherWright_internal( ...
%     D, s, mu, two_side_flag, iters, init_str)
function [freq_struct, absorption_struct, simulation_struct, q, ...
            num_absorptions, num_fixations, num_losses, ...
            num_absorptions_by_generation_vec, num_fixations_by_generation_vec, num_losses_by_generation_vec, ...
            count_vec] = simulate_forward_FisherWright_internal(D, s, mu, two_side_flag, iters, init_str)

if(~isfield(D, 'cond_on_polymorphic_flag'))
    D.cond_on_polymorphic_flag=1; % NEW!!! DEFAULT IS SIMULATE ONLY POLYMORPHIC ALLELES !!!
end
N_vec = demographic_parameters_to_n_vec(D, D.index);  % compute population size at each generation
N = N_vec(1); max_N = max(N_vec); num_generations = length(N_vec)-1;
if(~isfield(D, 'add_new_alleles') || isempty(D.add_new_alleles)) % default: add newly born alleles at each generation
    D.add_new_alleles = 1;
end
if(~isfield(D, 'compute_absorb')) % default: don't compute absorption time (can be heavy)
    %    if(nargout > 11) % compute absorption time and count vec. Can be heavy (?)
    %        D.compute_absorb = 1;
    %    else
    D.compute_absorb = 0;
    %    end
end
freq_struct.p_vec = cell(num_generations+1, 1); freq_struct.x_vec = freq_struct.p_vec; % het_vec = p_vec; % initilize distributions

rand_str = 'poisson'; % 'poisson'; % 'binomial'; % 'poisson'; % 'binomial'; % How to simulate each generation: poisson is much faster (approximation)
max_num_alleles = 20000; % maximum number of alleles to simulate (to save time!)
block_size = min(max_num_alleles, iters); % number of simulations to perform simultaniously
%final_q = zeros(iters, num_generations, 'single'); % fill this with alleles not absorbed
num_simulated_polymorphic_alleles_vec = zeros(num_generations+1, 1); % count how many iterations are left at each generation

% New: track only alleles that reached absorption (to avoid bias)
num_absorptions = 0; num_fixations = 0; num_losses = 0;
num_absorptions_by_generation_vec = zeros(num_generations,1);
num_fixations_by_generation_vec = num_absorptions_by_generation_vec;
num_losses_by_generation_vec = num_absorptions_by_generation_vec;
[absorption_struct.absorption_time_given_init_freq_vec, ...
    absorption_struct.fixation_time_given_init_freq_vec, ...
    absorption_struct.loss_time_given_init_freq_vec] = deal(zeros(2*max_N+1,1));


count_vec = zeros(2*max_N+1,1); % count how many alleles were at each allele frequency (from alleles absorbed)

num_alleles_simulated=0; total_polymorphic_generations = 0; ctr_alleles_blocks_simulated = 0;
mean_time_allele_polymorphic_at_equilibrium = absorption_time_by_selection(s, 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! time until absorbtion for a newly born allele
prob_site_polymorphic_at_equilibrium = (2*N*mu) * 2 * absorption_time_by_selection(s, 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start
total_het_at_each_generation_vec = zeros(num_generations, 1, 'single');
weights = [];

[generation_num_alleles, generation_weight] = num_alleles_to_simulate_per_generation_internal(N_vec, s, iters); % compute number of alleles to simulate (importance sampling)
generation_weight_expanded = RunLength(generation_weight, generation_num_alleles);

% TEMP: compute for small example also absorption probability at stationary
% M = FisherWright_ComputeMarkovMatrix(N_vec(1), s, 'exact', 1); % compute Markov matrix
% M(1,2)=mu; M(1,1)=1-mu; M(end,2)=mu; M(end,end)=1-mu; % add small mutation to make process ergodic

%[T_vec, T_mat] = MarkovChainAbsoptionTime(M, [1 2*N+1]); 
% pi_vec = markov_chain_stationary_dist(M);  % pi_poly_vec = [0 pi_vec(2:end-1)' ./ sum(pi_vec(2:end-1)) 0]'

q = zeros(iters, num_generations+1, 'single'); % NEW approach !! allocate in advance entire q array (should be faster) !!!
survive_prob = ones(size(q), 'single'); survive_prob_by_gen = zeros(num_generations, 1);  % ones(size(q(block_inds,:))); %[first_time_vec, last_time_vec] = deal(ones(1, iters)); % last_time_vec = ones(1, iters); % For each allele record the first and last polymorphic times
sim_ind_vec = zeros(iters, 1); % indicators saying if each coordinate was already simulated
while( (num_alleles_simulated < iters) && (num_simulated_polymorphic_alleles_vec(1) < max_num_alleles) ) % simulate blocks. Problem: No new alleles born here! (these can be simulated separatey?)
    block_inds = find(sim_ind_vec==0, min(generation_num_alleles(1), min(block_size, iters-sum(sim_ind_vec)))); sim_ind_vec(block_inds)=1; % set block indices
    %    weights = ones(block_size, num_generations, 'single');  % give later generations higher weights (why?)
    switch init_str % determine starting allele frequencies
        case 'equilibrium'
            q(block_inds,1) = single(round(2*N.* vec2column(allele_freq_spectrum_rnd(s, N, two_side_flag, length(block_inds))))); % sample allele frequency from equilibrium distribution
        case 'newly_born'
            q(block_inds,1) = ones(length(block_inds), 1, 'single'); % start with newly born alleles
    end
    num_simulated_polymorphic_alleles_vec(1) = num_simulated_polymorphic_alleles_vec(1)+ block_size;
    total_polymorphic_generations = total_polymorphic_generations + block_size; % + ITERS
    [unique_time_cum, rand_time_cum, arrange_time_cum] = deal(0);
    for j=1:num_generations % -1 % run vectorized. This is heaviest loop!
        if(mod(j,100)==0)
            fprintf('Sim. Wright-Fisher Gen. %ld out of %ld, #alleles=%ld\n', j, num_generations, length(block_inds));
            fprintf('      init-time=%.2f, rand-time=%.2f, arrange-time=%.2f\n', unique_time_cum, rand_time_cum, arrange_time_cum);
        end
        unique_time=cputime;
        if(isempty(weights))
            [U, C] = unique_with_counts(vec2row(q(block_inds,j))); % Compute histogram of counts
            total_het_at_each_generation_vec(j) = total_het_at_each_generation_vec(j) + ...
                2 .* sum(q(block_inds,j) ./ (2.*N_vec(j)) .* (1-q(block_inds,j) ./ (2.*N_vec(j)) ) ); % this indicates how much heterozygosity was absorbed at each time - should go down if we start at equilibrium!
        else  % here use weighted sums
            [U, C] = unique_with_counts(vec2row(q(block_inds,j)), [], weights); % Compute histogram of counts with weights !!!
            total_het_at_each_generation_vec(j) = total_het_at_each_generation_vec(j) + ...
                2 .* sum(vec2column(weights) .* q(block_inds,j) ./ (2.*N_vec(j)) .* (1-  q(block_inds,j) ./ (2.*N_vec(j)) ) ); % this indicates how much heterozygosity was absorbed at each time - should go down if we start at equilibrium!
        end
        % NEW !!! count by weights !!!
        [freq_struct.x_vec{j}, freq_struct.p_vec{j}] = union_with_counts(freq_struct.x_vec{j}, freq_struct.p_vec{j}, ...
            [0 U 2*N_vec(j)], [num_losses C num_fixations]);
        unique_time=cputime-unique_time; unique_time_cum = unique_time_cum+unique_time; rand_time=cputime; % get times
        
        %        new_block_inds = block_inds; % set indices for next generation
        if(D.add_new_alleles) % NEW! Add newly born alleles. Rate DOES NOT depend on current polymorphic probability !!! % Separate: simulate newly born alleles
            newly_born_inds = find(sim_ind_vec==0, min(generation_num_alleles(j+1), iters-sum(sim_ind_vec))); sim_ind_vec(newly_born_inds)=1; % set block indices;
        end
        [q([vec2row(block_inds) vec2row(newly_born_inds)], j+1), new_survive_prob, ~, num_simulated_polymorphic_alleles_vec(j+1)] = ... % q([vec2row(block_inds) vec2row(newly_born_inds)],j+1)
            simulate_WrightFisher_one_step(q(block_inds,j), N_vec, rand_str, s, j, ...
            D.cond_on_polymorphic_flag, D.add_new_alleles, num_simulated_polymorphic_alleles_vec(j+1), ...
            max_num_alleles, generation_num_alleles, generation_weight);
        survive_prob_by_gen(j) = mean(new_survive_prob(block_inds)); % survive_prob(block_inds,j));
        % Simulate one step (heavy part !!!). Always keep same size (block inds)
        if(D.add_new_alleles) % NEW! Add newly born alleles. Rate DOES NOT depend on current polymorphic probability !!! % Separate: simulate newly born alleles
            block_inds = [vec2row(block_inds) vec2row(newly_born_inds)];
        end
        survive_prob(block_inds,j+1) = survive_prob(block_inds,j) .* new_survive_prob; % multiply to get survive probs
        loss_inds = block_inds(find(q(block_inds,j+1) == 0)); fixation_inds = block_inds(find(q(block_inds,j+1) == 2*N_vec(j+1)));
        absorption_inds = vec2row(union(loss_inds, fixation_inds)); % reached fixation/extinsion and stop
        survived_inds = setdiff(block_inds, absorption_inds); % all alleles except absorbed ones !!!
        total_polymorphic_generations=total_polymorphic_generations+length(survived_inds);
        num_simulated_polymorphic_alleles_vec(j+1) = num_simulated_polymorphic_alleles_vec(j+1) + length(survived_inds);
        %last_time_vec(absorption_inds) = j; % set the time at which these indices were absorbed
        rand_time=cputime-rand_time; rand_time_cum = rand_time_cum+rand_time;  arrange_time=cputime;
        
        if(D.cond_on_polymorphic_flag) % temp: update weights
            weights = survive_prob(block_inds,j+1) .* generation_weight_expanded(block_inds); % NEW! weight alleles by their survival probabilities !!! (but what if we already had weights?)                       
        end
        
        if(D.compute_absorb) % compute absorption time and count vec. Can be heavy (?)
            if(~isempty(absorption_inds))
                for k=1:j % Alternative: loop on generations (not on indices of iterations)
                    [unique_inds, unique_counts] = unique_with_counts(  q(absorption_inds, k) );
                    unique_counts = unique_counts(unique_inds>0); unique_inds = unique_inds(unique_inds>0);
                    absorption_struct.absorption_time_given_init_freq_vec(unique_inds) = ...
                        absorption_struct.absorption_time_given_init_freq_vec(unique_inds) + (j-k+1) .* unique_counts;
                    count_vec(unique_inds) = ...
                        count_vec(unique_inds) + unique_counts;
                end
            end
        end
        % here we added q <- [q new_q], and then added new alleles !!!
        
        %        first_time_vec(absorption_inds) = j; % update for next time. No! we throw away the ones absorped !!
        num_absorptions_by_generation_vec(j) = num_absorptions_by_generation_vec(j) + length(absorption_inds); % count absorptions to see how much of the distirbution is kept
        num_fixations_by_generation_vec(j) = num_fixations_by_generation_vec(j) + length(fixation_inds);
        num_losses_by_generation_vec(j) = num_losses_by_generation_vec(j) + length(loss_inds);
        arrange_time=cputime-arrange_time;
        arrange_time_cum = arrange_time_cum+arrange_time;
        
        % TEMP DEBUG: 
%         p_init_vec = p_vec{1} ./ sum(p_vec{1})
%         p_next_vec = p_init_vec * M
%         p_abs = p_next_vec(1)+p_next_vec(end)
%         M * p_init_vec'
%         p_next_poly_vec = pi_poly_vec' * M
%         p_poly_abs = p_next_poly_vec(1)+p_next_poly_vec(end); 
%         M_empiric = zeros(2*N_vec(1)-1); 
%         for k=1:2*N_vec(1)-1
%             for m=1:2*N_vec(1)-1
%                 M_empiric(k,m) = sum((q(1:generation_num_alleles(1),1) == k) & (q(1:generation_num_alleles(1),2) == m)); 
%             end
%         end
%         M_empiric = M_empiric ./ repmat(sum(M_empiric,2), 1, 2*N_vec(1)-1); % normalize 
%         M_normalized = M(2:end-1,2:end-1); M_normalized = M_normalized ./ repmat(sum(M_normalized, 2), 1, 2*N_vec(1)-1); 
%         figure; imagesc(M_empiric-M_normalized); colorbar;
    end % loop on generations (heaviest loop!)
    fprintf('init-time=%f, rand-time=%f, arrange-time=%f\n', unique_time, rand_time_cum, arrange_time_cum)
    cur_num_alleles_survived = min(size(q, 1), iters - num_alleles_simulated);
    %     if(size(q,1) > iters  - num_alleles_simulated) % adjust last block
    %         num_simulated_polymorphic_alleles_vec = num_simulated_polymorphic_alleles_vec - size(q,1) ..; % adjust num alleles vec
    %     end
    
    %    final_q((num_alleles_simulated+1):(num_alleles_simulated+cur_num_alleles_survived),:) = ...
    %        q(1:cur_num_alleles_survived,:);
    num_alleles_simulated = num_alleles_simulated + cur_num_alleles_survived
    
    % compute again for last generation
    if(~exist('weights', 'var') || isempty(weights))
        [U, C] = unique_with_counts(vec2row(q(:, num_generations+1))); % Compute histogram of counts
    else % here use weighted sums
        [U, C] = unique_with_counts(vec2row(q(:, num_generations+1)), [], weights); % Compute histogram of counts with weights !!!
    end
    [freq_struct.x_vec{num_generations+1}, freq_struct.p_vec{num_generations+1}] = ...
        union_with_counts(freq_struct.x_vec{num_generations+1}, freq_struct.p_vec{num_generations+1}, ...
        [0 U 2*N_vec(num_generations)], [num_losses C num_fixations]);
    ctr_alleles_blocks_simulated = ctr_alleles_blocks_simulated + block_size;
end % while num_alleles <= iters

total_het_at_each_generation_vec = total_het_at_each_generation_vec ./ ...
    num_simulated_polymorphic_alleles_vec(1); % normalize by number of iterations
num_absorptions = sum(num_absorptions_by_generation_vec); % count absorptions, losses and fixations
num_fixations = sum(num_fixations_by_generation_vec);
num_losses = sum(num_losses_by_generation_vec);

for j=1:(num_generations+1) % sort and normalize distributions
    [freq_struct.x_vec{j}, sort_perm] = sort(freq_struct.x_vec{j});
    freq_struct.p_vec{j} = freq_struct.p_vec{j}(sort_perm) ./ num_simulated_polymorphic_alleles_vec(j); % new! normalize to sum to 1!!. (Before: (1) instead of (j):  normalize to incorporate prob. of alleles which started as polymorphic
end
%q=final_q; % Copy results - why?

switch init_str % Normalize
    case 'equilibrium'
        total_het_at_each_generation_vec = total_het_at_each_generation_vec .* prob_site_polymorphic_at_equilibrium; % normalize to include fraction of polymorphic sites at the beginning equilibrium
    case 'newly_born'
        total_het_at_each_generation_vec = total_het_at_each_generation_vec .* 2*N*mu; % normalize to include fraction of polymorphic sites born at each generation
        total_het_at_each_generation_vec = cumsum(total_het_at_each_generation_vec); % unite the contribution from each generation
end
prob_site_polymorphic_at_end = prob_site_polymorphic_at_equilibrium .* ... % esitmate Prob. poly, by ratio of # alleles at start and at end (noisy!)
    num_simulated_polymorphic_alleles_vec ./ num_simulated_polymorphic_alleles_vec(1);
absorption_struct.absorption_time_given_init_freq_vec = absorption_struct.absorption_time_given_init_freq_vec ./ count_vec; % normalize to get estimated absorption times
q = q ./ repmat(2.*N_vec(1:(end))', iters, 1); % transfer from counts to frequencies
frac_polymorphic_vec = 1-num_absorptions_by_generation_vec ./ iters;


% Now compute correction factor: this is the ratio between #of simulated
% alleles, and #alleles we would have expected for this simulation
L_correction_factor = 4*N*mu*mean_time_allele_polymorphic_at_equilibrium / num_alleles_simulated;

simulation_struct = var2struct(q, num_simulated_polymorphic_alleles_vec, weights, ...
    frac_polymorphic_vec, L_correction_factor);

% Simulate one step of Wright-Fisher process forward
% Input:
% old_q - vector of current allele frequencies
% N_vec - vector of population sizes
% rand_str - how to sample next generation
% s - selection coefficient
% j - index of current generation
% cond_on_polymorphic_flag - simulate only polymorphic alleles (default: 'off')
% add_new_alleles - add newly born alleles at each generation
% num_simulated_polymorphic_alleles - ???
% max_num_alleles - maximum number of alleles to simulate
% generation_num_alleles - number of new alleles to generate in each generation
% generation_weight - weight assigned for each generation (importance sampling)
%
% Output:
% new_q - allele frequencies in next generation
% survive_prob - survival probability when we force polymorphic alleles
% weights - give weight to each simulated allele
% num_simulated_polymorphic_alleles - total number of polymorphic alleles simulated
%
function [new_q, survive_prob, weights, num_simulated_polymorphic_alleles] = ...
    simulate_WrightFisher_one_step(old_q, N_vec, rand_str, s, j, ...
    cond_on_polymorphic_flag, add_new_alleles, num_simulated_polymorphic_alleles, ...
    max_num_alleles, generation_num_alleles, generation_weight) % Simulate one step (heavy part !!!))

if(~exist('cond_on_polymorphic_flag', 'var') || isempty(cond_on_polymorphic_flag)) % default: don't condition on polymorphic
    cond_on_polymorphic_flag=0;
end
new_expected_q = old_q .* ((1+s)./(1+s.*old_q./(2*N_vec(j)))) ./ (2*N_vec(j));  % new allele freq. of the deleterious alleles
switch rand_str % sample new generation
    case {'binomial', 'exact'}
        survive_prob =  1-binopdf(0, 2*N_vec(j+1), new_expected_q)-binopdf(2*N_vec(j+1), 2*N_vec(j+1), new_expected_q);
        mono_inds = 1:length(new_expected_q); new_q = zeros(size(old_q));
        while 1 % Simulate only polymorphic
            new_q(mono_inds) = binornd(2*N_vec(j+1), new_expected_q(mono_inds)); % ./ (2*N_vec(j+1)); % randomize next generation
            mono_inds = find((new_q == 0) | (new_q == 2*N_vec(j+1)));
            if(~cond_on_polymorphic_flag || isempty(mono_inds))
                break
            end
        end
    case {'poisson', 'approx', 'approximate'} % use poisson approximation (good for large N)
        % here change so that poisson/Gaussian criteria is the same as in numeric !!!
        small_inds = find(new_expected_q .* 2*N_vec(j+1) < 50); % low frequency alleles
        big_inds = setdiff(find((1-new_expected_q) .* 2*N_vec(j+1) < 50), small_inds); % for high frequency alleles poisson approximation isn't good enough and we simulate binomial distribution
        medium_inds = setdiff(1:length(new_expected_q), union(small_inds, big_inds));
        survive_prob = ones(size(new_expected_q));
        survive_prob(small_inds) = 1-poisspdf(0, double(2*N_vec(j+1) .* new_expected_q(small_inds))); % neglect fixation
        survive_prob(big_inds) = 1-poisspdf(2*N_vec(j+1), double(2*N_vec(j+1) .* new_expected_q(big_inds))); % neglect loss
        new_q = zeros(size(old_q));
        if(~isempty(medium_inds)) % for these alleles, Gaussian approximation: simulation does depend on N.
            survive_prob(medium_inds) = 1-normcdf(0, 2*N_vec(j+1) .* new_expected_q(medium_inds), ...
                sqrt(2*N_vec(j+1) .* new_expected_q(medium_inds) .* (1-new_expected_q(medium_inds))));  % neglect loss
            new_q(medium_inds) = round( normrnd( 2*N_vec(j+1) .* new_expected_q(medium_inds), ...
                sqrt(2*N_vec(j+1) .* new_expected_q(medium_inds) .* (1-new_expected_q(medium_inds))) ) ); % ./ (2*N_vec(j+1)); % randomize next generation
        end
        while 1 % update for small indices (close to 0) and large indices (close to 2N)
            new_q(small_inds) = poissrnd(double(2*N_vec(j+1) .* new_expected_q(small_inds))); % ./ (2*N_vec(j+1)); % randomize next generation
            new_q(big_inds) = 2*N_vec(j+1) - poissrnd(double(2*N_vec(j+1) .* (1-new_expected_q(big_inds))));
            if(~cond_on_polymorphic_flag || (isempty(small_inds) && isempty(big_inds))) % simulated all indices
                break
            end
            small_inds = small_inds(new_q(small_inds) == 0); % monomorphic !!
            big_inds = big_inds(new_q(big_inds) == 2*N_vec(j+1));  % these are monomorphic !!
            new_q = max(0, min(new_q, 2*N_vec(j+1))); % correct for Gaussian approximation giving negative results
            small_inds = union(small_inds, find(new_q == 0));
            big_inds = union(big_inds, find(new_q == 2*N_vec(j+1)));
        end
end % switch rand_str

if(~cond_on_polymorphic_flag)
    survive_prob = new_q>0; % determine who survived (0/1 probabilities)
    absorption_inds = find((new_q == 0) | (new_q == 2*N_vec(j+1)));
else
    absorption_inds = [];
end
survived_inds = setdiff(1:length(new_q), absorption_inds);

weights = []; % dummy

% New: add newly born alleles !!!
if(add_new_alleles) % NEW! Add newly born alleles. Rate DOES NOT depend on current polymorphic probability !!!
    new_q(absorption_inds) = 0; % set to zero ALL (why not just current generation j?)
    num_new_alleles = generation_num_alleles(j+1); % NEW!!! poissrnd( (block_size /(2*mean_time_allele_polymorphic_at_equilibrium)) * (N_vec(j+1)/N) ); % Proportional to mutation rate times # of chromosomes . Mutation rate is cancelled !
    num_simulated_polymorphic_alleles = num_simulated_polymorphic_alleles+num_new_alleles; % update BEFORE down-weighting!!! to keep weighted # !!!
    
    survive_prob(absorption_inds) = 0; survive_prob = [survive_prob' ones(1, num_new_alleles)]'; % add survival probs for new alleles  %-length(absorption_inds))]';
    new_q = [new_q' ones(1, num_new_alleles)]'; % add new alleles
else % don't add new alleles
    q = q(survived_inds,:); % take only indices that are left
end % if D.add_new_alleles


% Determine total number of alleles to simulate at each genertion and how to weight them
% Input:
% N_vec - population size at each generation
% s - selection coefficient
% iters - total number of alleles to simulate
%
% Output:
% generation_num_alleles - number of new alleles to simulate in each generation
% generation_weight - weight of alleles simulated in each generation
%
function [generation_num_alleles, generation_weight] = num_alleles_to_simulate_per_generation_internal(N_vec, s, iters)

%num_gen = length(N_vec); %num_new_alleles = poissrnd( (block_size /(2*mean_time_allele_polymorphic_at_equilibrium)) * (N_vec(j+1)/N) ); % Proportional to mutation rate times # of chromosomes . Mutation rate is cancelled !
if(N_vec(1)>100) % use approximation for large N 
    mean_time_allele_polymorphic_at_equilibrium = 1+2*absorption_time_by_selection(s, 1, N_vec(1), 1/(2*N_vec(1)), 0.999999999, 0);  % NEW! Factor of two here! time until absorbtion for a newly born allele    
else
    mean_time_allele_polymorphic_at_equilibrium = absorption_time_by_selection_numeric(s, N_vec(1));  % NEW! Factor of two here! time until absorbtion for a newly born allele
end

total_het_by_gen_vec = n_vec_to_het_per_generation(N_vec); 
% set relative proportions. First is old alleles (equilibrium).
generation_num_alleles = [1*mean_time_allele_polymorphic_at_equilibrium+1*2*0 vec2row(N_vec(2:end)) ./ N_vec(1)]'; % take only from gen. 2. Hack! Add 1 to equilibrium alleles (good for s=0)!!!


% TEMP! get exact absorption probability to determine ratio: 
% TEMP: compute for small example also absorption probability at stationary
% % % M = FisherWright_ComputeMarkovMatrix(N_vec(1), s, 'exact', 1); % compute Markov matrix
% % % mu=1; M(1,2)=mu; M(1,1)=1-mu; M(end,2)=mu; M(end,end)=1-mu; % add small mutation to make process ergodic
% % % 
% % % %[T_vec, T_mat] = MarkovChainAbsoptionTime(M, [1 2*N+1]); 
% % % pi_vec = get_stationary_dist(M);  pi_poly_vec = [0 pi_vec(2:end-1)' ./ sum(pi_vec(2:end-1)) 0]'


% New! take into account heterozygosity (importance sampling!) 
%generation_num_alleles = generation_num_alleles .* vec2row(total_het_by_gen_vec) ./ mean(total_het_by_gen_vec); 
% generation_weight = generation_weight ./ 

% Compute heterozygosity expected by site 

MIN_GEN = 1; % Make sure that each generation have enough alleles born at each generation
new_generation_num_alleles = max(MIN_GEN, ceil(iters * generation_num_alleles ./ ...
     sum(generation_num_alleles))); % set total number of alleles. Only iters?
% new_generation_num_alleles = max(MIN_GEN, ceil(iters * generation_num_alleles .* total_het_by_gen_vec ./ ...
%     sum(generation_num_alleles .* total_het_by_gen_vec))); % set total number of alleles. Only iters?
while(sum(new_generation_num_alleles) > iters)
    k = sum(new_generation_num_alleles)-iters; [~, sort_perm] = sort(new_generation_num_alleles, 'descend'); 
    big_inds = intersect(sort_perm(1:k), find(new_generation_num_alleles>1)); % Reduce 1 to make sum exactly equal to iters
    new_generation_num_alleles(big_inds) = new_generation_num_alleles(big_inds)-1;
end

generation_weight =  generation_num_alleles ./ new_generation_num_alleles; % give weights: #alleles which should be vs. #alleles simulated
generation_weight = generation_weight ./ mean(generation_weight); % normalize

poisson_flag=0; % set actual number of alleles with/without sampling
if(poisson_flag)
    generation_num_alleles = poissrnd(new_generation_num_alleles);
else
    generation_num_alleles = new_generation_num_alleles;
end

% Compute how much heterozygosity is retained on expectation in every generation 
% Input: 
% N_vec - population size at each generation 
% 
% Output: 
% total_het_by_gen_vec - expected heterozygosity remaining at end from each generation 
function total_het_by_gen_vec = n_vec_to_het_per_generation(N_vec)
s=0; % should later compute for general s 
total_het_by_gen_vec = cumprod(1 - 1./(2.*N_vec(end:-1:1))); total_het_by_gen_vec = total_het_by_gen_vec(end:-1:1); % heterozygosity per 2N alleles for alleles born at generation t
% First generation should have higher weight: it is at equilibrium and has both 
% (i) more alleles, and (ii) higher frequencies !! 
%mean_time_allele_polymorphic_at_equilibrium = absorption_time_by_selection(s, 1, N_vec(1), 1/(2*N_vec(1)), 0.999999999, 0);  % NEW! Factor of two here! time until absorbtion for a newly born allele
total_het_by_gen_vec(1) = total_het_by_gen_vec(1) * N_vec(1)/log(2*N_vec(1)); %taken care of only (ii): 


% set relative proportions. First is old alleles (equilibrium).
% generation_num_alleles = [2*mean_time_allele_polymorphic_at_equilibrium vec2row(N_vec(2:end)) ./ N_vec(1)]; % take only from gen. 2


%% New: (10/2/2019): we should also compute approxiamte solution for large population size (coarse-grained). 
%% Look at finite-differences methods for differential equations 
% Compute allele frequency distribution numerically.
% Method: for each generation propagate the Markov chain
% NOTE: Selection s should be NEGATIVE for deletirious alleles
% Todo: Change output to have structures
%
% Input:
% N - initial population size (why needed if present in N_vec)
% s - selection coefficient. Should be NEGATIVE for detelitirious alleles
% mu - mutation rate
% two_side_flag - derived/minor allele frequency
% num_generations - total number of generations to simulate
% N_vec - population size at each generation
% max_N - maximal population size (?)
% init_str - how to start (equilibrium or new allele)
% compute_matrix - flags saying if to comptue Markov Matrix M 
% num_bins - NEW! allow approximate SFS. Do not compute x_vec and p_vec for
% every count from 0 to 2*N, but divide to bins (logarithmically) 
%
% Output:
% x_vec - vector of num. carriers at each generation
% p_vec - vector of probabilities at each generation
% total_het_at_each_generation_vec - total heterozygosity at each generation
% num_absorptions - number of absorptions (old alleles which dissappeared)
% num_fixations - number of fixations (old alleles which fixate in the population)
% num_absorptions_by_generation_vec - number of alleles which were absorped per-generation
% count_vec - vector counting the .. ???
% absorption_time_given_init_freq_vec - time to absorptions for each allele which fixed
% fixation_time_given_init_freq_vec - time to fixation for each allele which fixed
% num_losses - number of losses (old alleles which are lost in the population)
% loss_time_given_init_freq_vec - at what generation each loss occurs
% frac_polymorphic_vec - fraction of polymorphic alleles at each generation
% M - return also transition matrix
%
% [x_vec, p_vec, total_het_at_each_generation_vec, ...
%     num_absorptions, num_fixations, num_absorptions_by_generation_vec, count_vec, ...
%     absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, ...
%     num_losses, loss_time_given_init_freq_vec, frac_polymorphic_vec, prob_site_polymorphic_at_end, M, ...
%     mu_vec_analytic, mu_vec_equilibrium, L_correction_factor] = ...
%     compute_numeric_forward_FisherWright_internal(s, mu, two_side_flag, N_vec, init_str, compute_matrix)
function  [freq_struct, absorption_struct, simulation_struct, ...
            num_absorptions, num_fixations, ...
            num_absorptions_by_generation_vec, count_vec, ...
            num_losses, M] = compute_numeric_forward_FisherWright_internal( ...
            s, mu, two_side_flag, N_vec, init_str, compute_matrix, num_bins)

num_absorptions = 0; num_fixations = 0; num_losses = 0; % temp
num_absorptions_by_generation_vec = []; count_vec = [];
N = N_vec(1); num_generations = length(N_vec)-1; max_N = max(N_vec);

[freq_struct.x_vec, freq_struct.p_vec, freq_struct.het_vec] = deal(cell(num_generations+1, 1));  % structures for keeping track of alleles
freq_struct.total_het_at_each_generation_vec = zeros(num_generations, 1, 'single');

if(~exist('num_bins', 'var') || isempty(num_bins))
    freq_struct.x_vec{1} = (0:2*N) ./ (2*N); % vector of allele frequencies (include monomorphic alleles)
    init_p_vec = allele_freq_spectrum_numeric(freq_struct.x_vec{1}, s, N, two_side_flag, 'linear') *N ; % TEMP! factor N.  get stationary distrubution at equilibrium (initial condition)
else % new! set x_vec as center of bins 
    freq_struct.x_vec{1} = [(0:50) linspace(50, 2*N, num_bins-50)] ./ (2*N); % set bins 
    bin_sizes = [1 diff(freq_struct.x_vec{1})*2*N]
    init_p_vec = exp( allele_freq_spectrum(freq_struct.x_vec{1}, s, N, two_side_flag, 'log') ) .* bin_sizes; % get stationary distrubution at equilibrium (initial condition)
    init_p_vec(1) = N; % correction for 0 frequency 
end
    
M=[]; % Markov chain
switch init_str
    case 'newly_born'
        init_p_vec(:) = 0; init_p_vec(1) = 1; % Set allele freq. at 0. Contribution will come from 1/2N
    case 'equilibrium'
        init_p_vec(end) = 0; % 1 allele frequency
end % switch init_str
init_p_vec = init_p_vec .* 2.* max(10^(-10), mu); %  .* 2*N;  % multiply by theta. (If mu=0: no mutations, then we've got a problem ... set some mu>0)
init_p_vec(1) = 1 - sum(init_p_vec(2:end-1)); % set first value to complete distribution to sum to one. Also when using init?

new_x_mean_vec = freq_struct.x_vec{1} .* (1+s); % first get the mean vector of new #offspring
%std_vec = sqrt(2*N .* new_x_mean_vec .* (1-new_x_mean_vec)); % binomial standard deviation
if(~exist('mu', 'var') || isempty(mu))
    mu = sum(vec2row(init_p_vec) .* binopdf(0, 2*N, freq_struct.x_vec) ./ (2*N)); % compute mutation rate per site (assuming T=1 target size)
end
freq_struct.p_vec{1} = vec2column(init_p_vec); het_vec{1} = 2.* freq_struct.p_vec{1} .* vec2column(freq_struct.x_vec{1} .* (1-freq_struct.x_vec{1}));

if(compute_matrix)
    p_vec_M = freq_struct.p_vec; p_vec_M{1} = p_vec_M{1}(1:end-1); % don't include alternative (fixed) allele
end
for j=1:num_generations % main loop on generations
    t_gen = cputime;
    [freq_struct.x_vec{j+1}, freq_struct.p_vec{j+1}, M, skip_run_frac, max_range] = ...
        compute_numeric_forward_FisherWright_one_step(freq_struct.x_vec{j}, freq_struct.p_vec{j}, s, mu, N_vec, j, compute_matrix);
    het_vec{j+1} = 2.* freq_struct.p_vec{j+1} .* vec2column(freq_struct.x_vec{j+1} .* (1-freq_struct.x_vec{j+1})); % update heterozygosity
    if(mod(j, ceil(50000/N)) == 0)
        fprintf('Comp. Markov-Chain, gen=%ld, skip=%.2f, range=%ld, gen-time=%.2f\n', j, skip_run_frac, max_range, cputime - t_gen)
    end
    %    generation_time = cputime - t_gen
end % loop on generations
% p_vec_equilibrium_analytic = p_vec{1};

% Compute summary statistics
freq_struct.total_het_at_each_generation_vec  = sum_cell(het_vec); % freq_struct.total_het_at_each_generation_vec(j+1) = sum(het_vec{j+1});
simulation_struct.frac_polymorphic_vec = zeros(num_generations,1);
for j=1:num_generations
    simulation_struct.frac_polymorphic_vec(j) = sum(freq_struct.p_vec{j}(2:end)); % This should change if p_vec isn't normalized !!     %    frac_polymorphic_vec(j) = integral_hist(freq_struct.x_vec{j}(2:end)', p_vec{j}(2:end)); % This should change if p_vec isn't normalized !!
    freq_struct.p_vec{j}(2:end) = freq_struct.p_vec{j}(2:end) ./ simulation_struct.frac_polymorphic_vec(j);
    freq_struct.p_vec{j}(1) = 0; % set to zero
end
if(compute_matrix)
    M2 = M; M2([1 end],:) = 0; ones_vec = ones(2*N,1); ones_vec(1) = 0; ones_vec(end) = 0;
    ttt=cputime; absorption_struct.absorption_time_given_init_freq_vec = linsolve(M2-eye(2*N), -ones_vec); inverse_time = cputime - ttt
else
    absorption_struct.absorption_time_given_init_freq_vec = zeros(2*N+1,1);
end
% prob_fixation = 1; % THIS IS WRONG. HOW TO COMPUTE THIS?
fixation_time_given_init_freq_vec = zeros(size(absorption_struct.absorption_time_given_init_freq_vec));
loss_time_given_init_freq_vec = zeros(size(absorption_struct.absorption_time_given_init_freq_vec)); % Time conditioned on fixation/loss. NEED TO FILL THESE

prob_site_polymorphic_at_end2 = zeros(num_generations,1); % NEW! compute also with recursive formula
prob_site_polymorphic_at_end2(1) = simulation_struct.frac_polymorphic_vec(1);
for j=1:num_generations+1  %change output format (convention)
    freq_struct.p_vec{j} = vec2row(freq_struct.p_vec{j});
    new_vec(j) = (1-prob_site_polymorphic_at_end2(j)) * 2*N_vec(j)*mu;
    old_vec(j) = prob_site_polymorphic_at_end2(j).^1 .* sum(freq_struct.p_vec{j}(2:end) .* (1 - freq_struct.x_vec{j}(2:end).^(2*N_vec(j)) - (1-freq_struct.x_vec{j}(2:end)).^(2*N_vec(j)) ) );
    prob_site_polymorphic_at_end2(j+1) = (1-prob_site_polymorphic_at_end2(j)) * 2*N_vec(j)*mu + ...
        prob_site_polymorphic_at_end2(j).^1 .* sum(freq_struct.p_vec{j}(2:end) .* (1 - freq_struct.x_vec{j}(2:end).^(2*N_vec(j)) - (1-freq_struct.x_vec{j}(2:end)).^(2*N_vec(j)) ) );
end
save('numeric_file', 'new_vec', 'old_vec');
for j=1:num_generations+1  %change output format (convention)
    freq_struct.x_vec{j} = round(freq_struct.x_vec{j} .* (2*N_vec(j))); % get number of individual alleles (not fracion)
end
prob_site_polymorphic_at_end = simulation_struct.frac_polymorphic_vec; % (end-1); % NEED TO FILL !!!

% NEW! compute EXACTLY the heterozygosity (separate to new and old).
if(~exist('num_moments', 'var') || isempty(num_moments))
    num_moments = 5; % how many moments to compute
end
mu_vec_equilibrium =  zeros(num_moments, 1); % need  a recurrence formula here !!
theta = 1;
for k=1:num_moments
    mu_vec_equilibrium(k) = absorption_time_by_selection(s, theta, N, 0, 1, -k-1);  % at equilibrium compute: int_{f=0}^1   f(1-f) \psi_0(f) df
end

mu_vec_analytic = zeros(num_moments, num_generations); % Comptue only the zero-th moment
prod_het_loss = cumprod(1 - 1./(2.*N_vec));
mu0_vec_old = mu_vec_equilibrium(1) .* prod_het_loss;
mu0_vec_new = zeros(num_generations+1, 1);
mu0_vec_new(1) = (1/2)/(2*N_vec(1));
for j=1:num_generations
    mu0_vec_new(j+1) = (1/2)/(2*N_vec(1)) + mu0_vec_new(j) * (1-1/(2*N_vec(j+1)));
end
mu_vec_analytic(1,:) = mu0_vec_old(1:num_generations)+mu0_vec_new(1:num_generations); % SHOULDNT GO DOWN !!!

mean_time_allele_polymorphic_at_equilibrium = absorption_time_by_selection(s, 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! time until absorbtion for a newly born allele
simulation_struct.L_correction_factor = 4*N*mu*mean_time_allele_polymorphic_at_equilibrium; %  / num_alleles_simulated;
simulation_struct.num_simulated_polymorphic_alleles_vec = 10000; 

% NEW! Function for advancing one generation in Markov chain computation
% Input:
% x_vec - frequency values
% p_vec - probability of alleles being at each frequency
% s - selection coefficient
% mu - mutation rate 
% N_vec - population sizes 
% j - generation
% compute_matrix - 1: compute Markov transition matrix. 0: No (default) 
% num_bins - new ! enable coarse grain resolution
%
% Output:
% new_x_vec - frequency values at next generation
% new_p_vec - probability of alleles being at each frequency at next generation
% M - Markov Transition matrix (optional)
% skip_running_fraction_of_inds - ???
% max_range - values with positive probability to compute  
function [new_x_vec, new_p_vec, M, skip_running_fraction_of_inds, max_range] = ...
    compute_numeric_forward_FisherWright_one_step(x_vec, p_vec, s, mu, N_vec, j, compute_matrix, num_bins)

if(~exist('compute_matrix', 'var') || isempty(compute_matrix))
    compute_matrix = 0; 
end
if(compute_matrix)
    M = zeros(2*N_vec(j)+1, 2*N_vec(j+1)+1);
else
    M = [];
end
num_sigmas = 5; % keep only these many st.d. (to speed-up computation)
normpdf_vec = normpdf(-6:10^(-5):6); % Prepare in advance Gaussian density to save time. Current resolution: 1 million
prob_compute_method = 'approximate'; % 'approximate'; % 'exact'; % 'approximate'; %'exact'; % s'exact'; % ''; % 'exact';
if (~exist('num_bins', 'var') || isempty(num_bins))
    new_x_vec = (0:2*N_vec(j+1)) ./ (2*N_vec(j+1)); % new: include also 0 and 1 freqs. to get absolute values
else
    new_x_vec = [(0:50) linspace(50, 2*N_vec(j+1), num_bins-50)] ./ (2*N_vec(j+1)); % set bins 
    bin_sizes = [1 diff(freq_struct.x_vec{1})*2*N]
%    init_p_vec = exp( allele_freq_spectrum(freq_struct.x_vec{1}, s, N, two_side_flag, 'log') ) .* bin_sizes; % get stationary distrubution at equilibrium (initial condition)
%    init_p_vec(1) = N; % correction for 0 frequency 
end

new_x_mean_vec = x_vec .* (1+s); % first get the mean vector for next generation
std_vec = sqrt(2*N_vec(j+1) .* new_x_mean_vec .* (1-new_x_mean_vec));
new_p_vec  = zeros(1, 2*N_vec(j+1)+1);
left_range_vec = max(2, round(new_x_mean_vec.*2*N_vec(j+1) - num_sigmas .* std_vec)); % NEW! start from 2, not 1 !!!! % calculate range to compute binomial distribution
right_range_vec = 1+min(2*N_vec(j+1), round(new_x_mean_vec.*2*N_vec(j+1) + num_sigmas .* std_vec));
non_negligile_x_inds = find(p_vec > 10^(-9)); % take only states with non-negligible probability
skip_running_fraction_of_inds = (2*N_vec(j)+1-length(non_negligile_x_inds)) / (2*N_vec(j)+1); % see how much time we've saved
max_range = max(right_range_vec - left_range_vec);
for k=vec2row(non_negligile_x_inds) % 1:2*N_vec(j)+1 % -1 % loop on current allele frequency. This is the heaviest part
    cur_range_vec = left_range_vec(k):right_range_vec(k); % set output range vec    
    switch prob_compute_method
        case 'exact'   % Exact binomial computation
            new_p_vec(cur_range_vec) = ...
                new_p_vec(cur_range_vec) + p_vec(k)  .* ...
                binopdf(cur_range_vec-1, 2*N_vec(j+1), new_x_mean_vec(k)); % why take cur_range_vec - 1 ???
            if(compute_matrix)
                M(k,:) = binopdf(0:2*N_vec(j+1), 2*N_vec(j+1), new_x_mean_vec(k)); % Why column here? % -1
            end
        case {'approx', 'approximate'} % Use Poisson approximation
            if(new_x_mean_vec(k) * 2*N_vec(j+1) < 50)  % use poisson approximation  to binomial on left tail
                new_p_vec(cur_range_vec) = ...
                    new_p_vec(cur_range_vec) + p_vec(k)  .* ...
                    poisspdf(cur_range_vec-1, 2*N_vec(j+1) * new_x_mean_vec(k));  % Replace with poisson approximation for larger values!!!
                if(compute_matrix)
                    M(k,:) = poisspdf(0:2*N_vec(j+1), 2*N_vec(j+1) * new_x_mean_vec(k)); % Poisson approximation
                end
            else % for larger values use poisson or Gaussian approximations
                if((1-new_x_mean_vec(k)) .* 2*N_vec(j+1) < 50) %  use poisson approximation to binomial on right tail,
                    new_p_vec(cur_range_vec) = ...
                        new_p_vec(cur_range_vec) + p_vec(k)  .* ...
                        poisspdf(2*N_vec(j+1) - (cur_range_vec-1), 2*N_vec(j+1) * (1-new_x_mean_vec(k)));  % This part is slowest. We replaced with poisson !!!
                    if(compute_matrix)
                        M(k,:) = poisspdf(2*N_vec(j+1):-1:0, 2*N_vec(j+1) * (1-new_x_mean_vec(k))); % Poisson approximation
                    end
                else % use Gaussian approximation to binomial in the middle of distribution
                    new_p_vec(cur_range_vec) = ...
                        new_p_vec(cur_range_vec) + p_vec(k)  .* ...
                        normpdf_vec( round(( 6+(cur_range_vec-1 - 2*N_vec(j+1) * new_x_mean_vec(k)) ./ std_vec(k) ) .* 10^5) ) ./ std_vec(k);
                    %                 normpdf(cur_range_vec-1, 2*N_vec(j+1) * new_x_mean_vec(k), ...
                    %                 sqrt(2*N_vec(j+1) * new_x_mean_vec(k) * (1-new_x_mean_vec(k))) );  % This part is slowest. We replace with poisson !!!
                    if(compute_matrix)
                        M(k,:) = normpdf(0:2*N_vec(j+1), 2*N_vec(j+1)*new_x_mean_vec(k), ...
                            sqrt(2*N_vec(j+1)*new_x_mean_vec(k)* (1-new_x_mean_vec(k))) ); % Why column here? % -1
                    end
                end
            end % if: what probability propagating approximation to use
            if (exist('num_bins', 'var') && ~isempty(num_bins))
                new_p_vec(cur_range_vec) = new_p_vec(cur_range_vec) .* bin_sizes(cur_range_vec)
            end
    end % switch approximate
end % loop on k
new_p_vec(1) = 1-sum(new_p_vec(2:end-1)); new_p_vec(end) = 0; % combine zero and one frequencies
new_p_vec(2) = new_p_vec(2) + (1*new_p_vec(1)+0) .* 2*N_vec(j+1)*mu*1; % add mutations. Try factor 2 !!
new_p_vec(1) = new_p_vec(1) .* (1-2*N_vec(j+1)*mu*1); % reduce monomorphic zero alleles due to mutations . Try factor 2 !!!
new_p_vec = vec2column(new_p_vec); %    p_vec{j+1} = vec2column(p_vec{j+1} ./ sum(p_vec{j+1})); % normalize to sum to one. Why not normalize at end???
if(compute_matrix)
    M(:,1) = M(:,1) + M(:, 2*N_vec(j+1)+1);
    M = M(1:(2*N_vec(j)), 1:(2*N_vec(j+1))); % identify ancestral and derived !!!
    M(:,2) = M(:,2) + M(:,1) .* 2*N_vec(j+1)*mu; % correct for mutation
    M(:,1) = M(:,1) .* (1-2*N_vec(j+1)*mu);
end



% Compute FisherWright SFS analytically (how? look at formulas for moments)
% NOT WORKING YET FOR NON-EQUILIBRIUM !!!
%
% Input:
% N - population size
% s - selection coefficient. NOTE: Should be NEGATIVE for deletirious alleles
% mu - mutation rate
% two_side_flag - derived/minor allele frequency
% num_generations - # generations with expansion
% N_vec - population size at each generation
% max_N - maximum possible populatin size
% init_str - how to start (equilibrium vs. new allele)
%
% Output:
% lots of things
%
% mu_vec_analytic - moments
% mu_vec_equilibrium - moments for equilibrium
%
function [freq_struct, absorption_struct, simulation_struct, ...
            num_absorptions ,num_fixations, ...
            num_absorptions_by_generation_vec, count_vec, num_losses] = ...
            compute_analytic_forward_FisherWright_internal( ...
            s, mu, two_side_flag, num_generations, N_vec, init_str, num_moments)
use_moments = 1;
use_gegenbauer = 0;
N = N_vec(1); max_N = max(N_vec);

[total_het_at_each_generation_vec, num_absorptions ,num_fixations, ...
    num_absorptions_by_generation_vec, count_vec ...
    absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, num_losses, ...
    loss_time_given_init_freq_vec, frac_polymorphic_vec] = deal([]);
x_vec = cell(num_generations, 1); p_vec = x_vec;
freq_struct.prob_site_polymorphic_at_equilibrium = (2*N*mu) * 2 * ...
    absorption_time_by_selection(s, 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start

if(use_gegenbauer) % this works only for constant population size
    max_deg = 20; % how many polynomials to sum
    coeff_vec = 2 .* (1:max_deg-1) .* (1-(1-2.*p)).^2 ./ ((1:max_deg).*(2:max_deg+1));
    gegenbauer_vec = zeros(max_deg, 1); % These work only for constant population size !
    for i=1:max_deg
        gegenbauer_vec(i) = mfun('G', i-1, 3/2, (1-2.*x_vec{1})) .* mfun('G', i-1, 3/2, (1-2.*p)); %   Gegenbauer polynomials
    end
    time_vec = exp(-(1:num_generations) ./ (4*N));
    for j=1:num_generations
        p_vec{j} = sum( time_vec(j).^ ((1:max_deg).*(2:max_deg+1)) .* coeff_vec .* gegenbauer_vec);
    end
end % use gegenbauer

if(use_moments) % this works only for s=0 !!!!
    freq_struct.prob_site_polymorphic_at_end = zeros(num_generations, 1);
    freq_struct.prob_site_polymorphic_at_end(1) = freq_struct.prob_site_polymorphic_at_equilibrium; % initialize
    %    prob_site_polymorphic_at_end3 = prob_site_polymorphic_at_end;
    if(~exist('num_moments', 'var') || isempty(num_moments))
        num_moments = 5; % how many moments to compute
    end
    [freq_struct.mu_vec_analytic, freq_struct.mu_vec_equilibrium] = ...
        FisherWright_Compute_SFS_Moments(N_vec, 0, num_moments, init_str); % compute moments with Formulas from Ewens
    for j=1:num_generations
        run_j = j
        x_vec{j} = (1:(2*N_vec(j)-1)) ./ (2*N_vec(j));
        % Estimate density using the max-entropy method:
        if(j==1)
            lambda0 = [];
        else % 'warm' start: take initial condition from previous generation
            lambda0 = lambda_max_ent;
        end
        [lambda_max_ent, g_het_max_ent, entr_max_ent] = ...  % Fit density using moments % , lambda0)
            me_dens2(mu_vec_analytic(2:end,j) ./ mu_vec_analytic(1,j), ...
            x_vec{j}, lambda0, 0); % normalize by 0th moment. Don't plot anything
        
        %         f_max_ent = zeros(1, 2*N_vec(j)-1);
        %         for k=1:num_moments
        %             f_max_ent = f_max_ent + lambda_max_ent(k) .* x_vec{j} .^ (k-1);
        %         end
        %         f_max_ent = exp(-f_max_ent) ./ (x_vec{j} .* (1-x_vec{j})); % Compute f. How to normalize?
        f_max_ent = mu_vec_analytic(1,j) .* g_het_max_ent' ./ (x_vec{j} .* (1-x_vec{j}));
        
        % New Normalization !!!
        %        prob_site_polymorphic_at_end2(j) = (N_vec(1)/N_vec(j)) * sum(f_max_ent * 2*mu)./2;
        
        
        
        
        
        % Yet another one! based on exponential distribution
        %         prob_site_polymorphic_at_end3(j+1) = (1-prob_site_polymorphic_at_end3(j)) * 2*N_vec(j)*mu;
        %         p_absorb3=1;
        %         for k=1:num_moments
        %             p_absorb3=p_absorb3-sum( f_max_ent .* ( (-2*N_vec(j+1).*x_vec{j}).^k ./ factorial(k)) ) / sum(f_max_ent);
        %             prob_site_polymorphic_at_end3(j+1) = prob_site_polymorphic_at_end3(j+1) - ...
        %                 prob_site_polymorphic_at_end3(j) * ...
        %                 sum( f_max_ent .* ( (-2*N_vec(j+1).*x_vec{j}).^k ./ factorial(k)) ) / sum(f_max_ent); %  UPDATE PROB. OF BEING POLYMORPHIC
        %         end
        
        
        f_max_ent = f_max_ent ./ sum(f_max_ent);
        %        p_absorb = sum ( f_max_ent .* (x_vec{j}.^(2*N_vec(j)) + (1-x_vec{j}).^(2*N_vec(j))) );
        %         if(j < num_generations)
        %             new_vec(j) = prob_site_polymorphic_at_end(j) * (1-p_absorb);
        %             old_vec(j) = (1-prob_site_polymorphic_at_end(j)) * 2*N_vec(j)*mu;
        %             prob_site_polymorphic_at_end(j+1) = prob_site_polymorphic_at_end(j) * (1-p_absorb) + ...
        %                 (1-prob_site_polymorphic_at_end(j)) * 2*N_vec(j)*mu; %  UPDATE PROB. OF BEING POLYMORPHIC
        %         end
        
        %        prob_site_polymorphic_at_end2(j+1) = (1-prob_site_polymorphic_at_end2(j)) * 2*N_vec(j)*mu + ...
        %            prob_site_polymorphic_at_end2(j).^1 .* sum(p_vec{j}(2:end) .* (1 - x_vec{j}(2:end).^(2*N_vec(j)) - (1-x_vec{j}(2:end)).^(2*N_vec(j)) ) );
        
        p_vec{j} = [0 f_max_ent 0]; % get p. Normalize and add zeros at boundaries
        
        freq_struct.prob_site_polymorphic_at_end(j) = 4*N_vec(1)*mu*mu_vec_analytic(1,j) / moment_hist(f_max_ent, ...
            f_max_ent .* vec2row(x_vec{j} .* (1-x_vec{j})), 0, 0, 0);
        x_vec{j} = [0 round(x_vec{j}*2*N_vec(j)) 2*N_vec(j)]; % get x in integers
        
        
    end % loop on number of generations
    %    prob_site_polymorphic_at_end = prob_site_polymorphic_at_end .* (N_vec(1)./N_vec(1:num_generations));
    %    prob_site_polymorphic_at_end = prob_site_polymorphic_at_end2';
end % use moments. Function doesn't use mutation rate mu !!

%save('moments_file', 'new_vec', 'old_vec');


% Temp: Plot individual trajectories of simulations
%
% Input:
% N - population size
% mu - mutation rate
% s - selection coefficient
% q - individual allele-frequency simulations
% p_vec - probability at each allele frequency
% p_vec_equilibrium_analytic - theoretical solution for population at equilibrium
%
function plot_expansion_internal(N, mu, s, q, x_vec, p_vec, p_vec_equilibrium_analytic)

params_str = ['N=' num2str(N) ', \mu=' num2str(mu) ', s=' num2str(s)];
figure; hold on; plot(q(1:20,:)'); xlabel('Time (generations)'); ylabel('Freq.');
title(['Trajectorues of Fisher-Wright model with selection and mutations. ' params_str]);
figure;  imagesc(q); colorbar;
title(['Trajectorues of Fisher-Wright model with selection and mutations. ' params_str]);
xlabel('Time (generations)'); ylabel('Iteration');

% Need to set x_vec:

figure; hold on;
plot(x_vec(2:(end-1)), p_vec);
plot(x_vec(2:(end-1)), p_vec_equilibrium_analytic, 'r');
xlabel('derived allele freq.'); ylabel('density');
legend('empirical', 'theoretical');
title(['Density of alleles for each allele frequency. ' params_str ...
    ' (' num2str(num_generations) ' generations, ' num2str(iters) ...
    ' iterations, ' num2str(num_effective_iters) ' effective iters)']);


% Compute how much old vs. new heterozygosity exists in sample
%
% Input:
% N - initial population size
% mu - mutation rate
% s - selection coefficient
% num_generations - number of generations for expansion
% expansion_factor - expansion rate at each generation
%
% Output:
% H_old - old heterozygosity
% H_new - new heterozygosity
%
%function [H_old, H_new] = het_expan(N, mu, s, num_generations, expansion_factor)
%
%N_vec = N .* expansion_factor.^(0:num_generations);
%prod_vec = cumprod( 1 - 1 ./ (2.*N_vec) );
%H_old = prod_vec(end); % prod ( 1 - 1 ./ (2.*N_vec) ); % assuming total heterozgosity is 1
%H_new = mu .* sum(prod_vec); % H_new_per_gen = mu;
%
