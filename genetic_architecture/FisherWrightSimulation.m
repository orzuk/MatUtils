% Simulate a generalization of Fisher-Wright model with selection and mutations.
% We allow expansion, and keep track of both old and new distribution of allele frequencies
%
% Input:
% region_flag - NEW! flag saying if to simulate a region with fixed size
%               (and then simulate iters such regions!) or simulate a fixed number of alleles (default)
% D - NEW! Demographic model. Contains population size at each generation !!!
% N - effective population size (either initial number or vector of population sizes at each generation)
% mu - mutation rate (per site/region per generation)
% s - selection coefficient. NOTE: Should be NEGATIVE for deletirious alleles
% num_generations - how many generations to run after we establish equilibrium
% expansion_factor - by how much does the population grow in each generation
% init_str - initialization of allele frequencies (equilibrium or newly-born alleles)
% iters - how many times to perform simulation
% compute_mode - how to compute (simulation/markov-chain numeric computation/analytic)
% num_bins - # of bins in returend histograms (when N is large we want a coarse-grained version of allele frequencies)
%
% Output:
% freq_struct - structure with output on allele frequencies. Fields:
%               x_vec - vector of allele frequencies for each generation
%               p_vec - vector of their probabilities (not counts: how many alleles present) for each generation.
%               p_vec_equlibrium_analytic - vector of allele frequencies computed analytically
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

absorption_time_given_init_freq_vec = []; q = []; % count_vec = [];

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
if(~isfield('index', D)) % assume only one demography
    D.index = 1;
end
if(~isfield('add_new_alleles', D)) % generate new alleles in each generation
    D.add_new_alleles = 1;
end


% Extract all infromation from Demographic model
N_vec = demographic_parameters_to_n_vec(D, D.index);  % compute population size at each generation
N = N_vec(1); num_generations = length(N_vec)-1; max_N = max(N_vec);
prob_site_polymorphic_at_equilibrium = (2*N*mu) * 2 * absorption_time_by_selection(abs(s), 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start

%x_vec = (0:2*N) ./ (2*N); % vector of allele frequencies (include monomorphic alleles) at start
%heterozygosity_per_site = 4*N*mu * absorption_time_by_selection(abs(s), 1, N, 0, 1, 'var');  % heterozygosity per site
%p_vec = cell(num_generations+1, 1); %total_het_at_each_generation_vec = zeros(num_generations,1);
switch compute_mode
    case 'simulation' % should be a sub-routine
        [q, weights, x_vec, p_vec, total_het_at_each_generation_vec, ...
            num_absorptions, num_fixations, num_losses, ...
            num_absorptions_by_generation_vec, num_fixations_by_generation_vec, num_losses_by_generation_vec, ...
            absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, loss_time_given_init_freq_vec, ...
            num_simulated_polymorphic_alleles_vec, count_vec, prob_site_polymorphic_at_end] = ...
            simulate_expansion_internal( ...
            N, s, mu, two_side_flag, iters, num_generations, N_vec, max_N, init_str, D.add_new_alleles);
        q = q ./ repmat(2.*N_vec(1:(end-1))', iters, 1); % transfer from counts to frequencies
        frac_polymorphic_vec = 1-num_absorptions_by_generation_vec ./ iters;
        num_effective_iters = sum(p_vec{end}(2:end-1)) % total number of single-generation single-allele steps performed (?)
        
    case 'numeric' % here compute everything by matrix multiplications
        [x_vec, p_vec, total_het_at_each_generation_vec, num_absorptions, num_fixations, ...
            num_absorptions_by_generation_vec, count_vec, ...
            absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, num_losses, ...
            loss_time_given_init_freq_vec, frac_polymorphic_vec, prob_site_polymorphic_at_end, M] = ...
            compute_numeric_expansion_internal( ...
            s, mu, two_side_flag, N_vec, init_str, compute_matrix);
        
    case {'analytic', 'moment', 'moments'} % compute solution based on Jacobi polynomials (See e.g. Kryukov et al. PNAS 2009).
        num_moments = min(floor(N_vec(1)/3), 5);
        [x_vec, p_vec, total_het_at_each_generation_vec, num_absorptions ,num_fixations, ...
            num_absorptions_by_generation_vec, count_vec ...
            absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, num_losses, ...
            loss_time_given_init_freq_vec, frac_polymorphic_vec, prob_site_polymorphic_at_end] = ...
            compute_analytic_expansion_internal( ...
            s, mu, two_side_flag, num_generations, N_vec, init_str, num_moments);
        
end % switch compute mode

% Fill additional statistics
het_vec = p_vec; % initilize distirbutions
for j=1:num_generations % heterozygosity vector
    het_vec{j} = 2 .* p_vec{j} .* vec2row(x_vec{j}./(2*N_vec(j)) .* (1-x_vec{j}./(2*N_vec(j))));
end
final_x_vec = (1:2*N_vec(num_generations)-1) ./ (2*N_vec(num_generations)); % set new coordinates
p_vec_equlibrium_analytic = exp( allele_freq_spectrum([0 final_x_vec 1], s, N, two_side_flag, 'log') ); % compute analytic approxiamtion (valid only for constant population size)
p_vec_equlibrium_analytic = normalize_hist(final_x_vec, p_vec_equlibrium_analytic(2:end-1)); % normalized



all_new_p_vec = []; all_new_x_vec = [];
switch init_str % unite distributions into one
    case 'newly_born' % we start at newly born alleles
        for j=1:num_generations
            [all_new_x_vec, all_new_p_vec] = union_with_counts(all_new_x_vec, all_new_p_vec, ...
                x_vec{j}, p_vec{j} .* 2*mu * N_vec(num_generations-j+1)); % compute weighted sum
        end
        [all_new_x_vec, sort_perm] = sort(all_new_x_vec);
        all_new_p_vec = all_new_p_vec(sort_perm);
        all_new_het_vec = 2.* all_new_p_vec .* ...
            all_new_x_vec./(2*N_vec(num_generations)) .* (1-all_new_x_vec/(2*N_vec(num_generations)));
    otherwise % we start at equilibrium
        all_new_het_vec = [];
        all_old_x_vec = x_vec{num_generations} .* (2*N_vec(num_generations));
        all_old_p_vec = p_vec{num_generations};
        all_old_het_vec = het_vec{num_generations};
end
prob_fixation = num_fixations / num_absorptions;
frac_old_alleles_survived_vec = [1 vec2row(cumprod(frac_polymorphic_vec))]; % take total # of alleles left
frac_het_vec = 99999; % this is the fraction of heterozygosity that is retained at each generation (TEMP!!!)
frac_het_kept_vec = [1 vec2row(cumprod(frac_het_vec))]; % take total # of alleles left
if(plot_flag) % Plot individual trajectories
    plot_expansion_internal(N, mu, s, q, x_vec, p_vec, p_vec_equlibrium_analytic);
end % if plot

% Build output structures
freq_struct = var2struct(x_vec, p_vec, p_vec_equlibrium_analytic, ...
    final_x_vec, het_vec, total_het_at_each_generation_vec, ...
    all_new_x_vec, all_new_p_vec, all_new_het_vec, compute_mode, ...
    prob_site_polymorphic_at_equilibrium, prob_site_polymorphic_at_end);
absorption_struct = var2struct(absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, loss_time_given_init_freq_vec, ...
    frac_polymorphic_vec, prob_fixation, frac_old_alleles_survived_vec, frac_het_kept_vec);
switch compute_mode % add simulations details
    case 'simulation'
         simulation_struct = var2struct(q, num_simulated_polymorphic_alleles_vec, ...
            frac_polymorphic_vec, prob_fixation, frac_het_kept_vec, ...
            num_simulated_polymorphic_alleles_vec); % add many things here
    otherwise % do not hold simulation information 
        simulation_struct = [];
end
simulation_time = cputime - simulation_time











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Internal functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute how much old and new heterozygosity exists
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
function [H_old, H_new] = het_expan(N, mu, s, num_generations, expansion_factor)

N_vec = N .* expansion_factor.^(0:num_generations);

prod_vec = cumprod( 1 - 1 ./ (2.*N_vec) );
H_old = prod_vec(end); % prod ( 1 - 1 ./ (2.*N_vec) ); % assuming total heterozgosity is 1
H_new = mu .* sum(prod_vec); % H_new_per_gen = mu;





%%%%%%%%%%%%%%%% Internal function for simulation
% Internal function for simulation of alleles to get frequency distribution.
% Simulate until we've got enough 'good alleles'
%
% Input:
% N - initial population size (number of INDIVIDUALS!)
% s - selection coefficient. Should be NEGATIVE for detelitirious alleles
% mu - mutation rate: per REGION per individual per generation
% two_side_flag - derived/minor allele frequency
% iters - how many alleles to output
% num_generations - total number of generations to simulate
% N_vec - population size at each generation
% max_N - maximal population size (?)
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
function [q, weights, x_vec, p_vec, total_het_at_each_generation_vec, ...
    num_absorptions, num_fixations, num_losses, ...
    num_absorptions_by_generation_vec, num_fixations_by_generation_vec, num_losses_by_generation_vec, ...
    absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, loss_time_given_init_freq_vec, ...
    num_simulated_polymorphic_alleles_vec, count_vec, prob_site_polymorphic_at_end] = ...
    simulate_expansion_internal( ...
    N, s, mu, two_side_flag, iters, num_generations, N_vec, max_N, init_str, add_new_alleles)


if(~exist('add_new_alleles', 'var') || isempty(add_new_alleles)) % default: add newly born alleles at each generation
    add_new_alleles = 1;
end

p_vec = cell(num_generations+1, 1); x_vec = p_vec; % het_vec = p_vec; % initilize distributions

rand_str = 'poisson'; % 'binomial'; % 'poisson'; % 'binomial'; % How to simulate each generation: poisson is much faster (approximation)
block_size = min(1000, iters); % number of simulations to perform simultaniously
final_q = zeros(iters, num_generations, 'single'); % fill this with alleles not absorbed
num_simulated_polymorphic_alleles_vec = zeros(num_generations, 1); % count how many iterations are left at each generation

% New: track only alleles that reached absorption (to avoid bias) 
num_absorptions = 0; num_fixations = 0; num_losses = 0;
num_absorptions_by_generation_vec = zeros(num_generations,1); 
num_fixations_by_generation_vec = num_absorptions_by_generation_vec; 
num_losses_by_generation_vec = num_absorptions_by_generation_vec;
absorption_time_given_init_freq_vec = zeros(2*max_N+1,1); fixation_time_given_init_freq_vec = zeros(2*max_N+1,1); loss_time_given_init_freq_vec = zeros(2*max_N+1,1);
count_vec = zeros(2*max_N+1,1); % count how many alleles were at each allele frequency (from alleles absorbed)

num_alleles_simulated=0; total_polymorphic_generations = 0; ctr_alleles_blocks_simulated = 0

mean_time_allele_polymorphic_at_equilibrium = absorption_time_by_selection(abs(s), 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! time until absorbtion for a newly born allele
prob_site_polymorphic_at_equilibrium = (2*N*mu) * 2 * absorption_time_by_selection(abs(s), 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start
total_het_at_each_generation_vec = zeros(num_generations, 1, 'single');

while( (num_alleles_simulated < iters) && (num_simulated_polymorphic_alleles_vec(1) < 20000) ) % simulate blocks. Problem: No new alleles born here! (these can be simulated separatey?)    
    q = zeros(block_size, num_generations, 'single');  % matrix of derived allele frequencies at each generation
    weights = ones(block_size, num_generations, 'single');  % give later generations higher weights (why?)
    switch init_str % determine starting allele frequencies
        case 'equilibrium'
            q(:,1) = round(2*N.* allele_freq_spectrum_rnd(s, N, two_side_flag, block_size)); % sample allele frequency from equilibrium distribution
        case 'newly_born'
            q(:,1) = 1; % start with newly born alleles
    end
    
    first_time_vec = ones(1, iters); last_time_vec = ones(1, iters); % For each allele record the first and last polymorphic times
    num_simulated_polymorphic_alleles_vec(1) = num_simulated_polymorphic_alleles_vec(1)+ block_size;
    total_polymorphic_generations = total_polymorphic_generations + iters;
    for j=1:num_generations-1 % run vectorized
        if(mod(j,100)==0)
            run_generation = j
        end
        %        cur_mu = mu .* N_vec(j+1) / N; % determine how many new mutations arise at each generation
        [U, C] = unique_with_counts(vec2row(q(:,j))); % Compute histogram of counts
        [x_vec{j}, p_vec{j}] = union_with_counts(x_vec{j}, p_vec{j}, [0 U 2*N_vec(j)], [num_losses C num_fixations]);
        new_q = q(:,j) .* ((1+s)./(1+s.*q(:,j)./(2*N_vec(j)))) ./ (2*N_vec(j));  % new allele freq. of the deleterious alleles
        total_het_at_each_generation_vec(j) = total_het_at_each_generation_vec(j) + ...
            2 .* sum(q(:,j) ./ (2.*N_vec(j)) .* (1-  q(:,j) ./ (2.*N_vec(j)) ) ); % this indicates how much heterozygosity was absorbed at each time - should go down if we start at equilibrium!
        
        switch rand_str % sample new generation
            case 'binomial'
                q(:,j+1) = binornd(2*N_vec(j+1), new_q); % ./ (2*N_vec(j+1)); % randomize next generation
            case 'poisson' % use poisson approximation (good for large N)
                big_inds = find(new_q>0.2); % for high frequency alleles poisson approximation isn't good enough and we simulate binomial distribution
                q(:,j+1) = poissrnd(double(2*N_vec(j+1) .* new_q)); % ./ (2*N_vec(j+1)); % randomize next generation
                if(~isempty(big_inds)) % for these alleles, simulation does depend on N.
                    q(big_inds,j+1) = normrnd( 2*N_vec(j+1) .* new_q(big_inds), ...
                        sqrt(2*N_vec(j+1) .* new_q(big_inds) .* (1-new_q(big_inds))) ); % ./ (2*N_vec(j+1)); % randomize next generation
                    q(big_inds,j+1) = max(0, round(q(big_inds,j+1)));
                end
                q(:,j+1) = min(q(:,j+1), 2*N_vec(j+1));
                if(mod(j,100)==0)
                    if(~isempty(big_inds))
                        frac_full_binomial_simulation = length(big_inds) / iters
                    end
                end
        end % switch rand_str
        loss_inds = find(q(:,j+1) == 0); fixation_inds = find(q(:,j+1) == 2*N_vec(j+1));
        absorption_inds = union(loss_inds, fixation_inds); % reached fixation/extinsion and stop
        survived_inds = setdiff(1:size(q,1), absorption_inds);
        total_polymorphic_generations=total_polymorphic_generations+length(survived_inds);
        num_simulated_polymorphic_alleles_vec(j+1) = num_simulated_polymorphic_alleles_vec(j+1) + length(survived_inds);
        last_time_vec(absorption_inds) = j; % set the time at which these indices were absorbed
        
        if(nargout > 11) % compute absorption time and count vec
            if(~isempty(absorption_inds))
                for k=1:j % Alternative: loop on generations (not on indices of iterations)
                    %                    cur_gen = j
                    %                    record_gen = k
                    [unique_inds, unique_counts] = unique_with_counts(  q(absorption_inds, k) );
                    unique_counts = unique_counts(unique_inds>0); unique_inds = unique_inds(unique_inds>0);
                    absorption_time_given_init_freq_vec(unique_inds) = ...
                        absorption_time_given_init_freq_vec(unique_inds) + (j-k+1) .* unique_counts;
                    count_vec(unique_inds) = ...
                        count_vec(unique_inds) + unique_counts;
                end
            end
        end
        %        first_time_vec(absorption_inds) = j; % update for next time. No! we throw away the ones absorped !!
        num_absorptions_by_generation_vec(j) = num_absorptions_by_generation_vec(j) + length(absorption_inds); % count absorptions to see how much of the distirbution is kept
        num_fixations_by_generation_vec(j) = num_fixations_by_generation_vec(j) + length(fixation_inds);
        num_losses_by_generation_vec(j) = num_losses_by_generation_vec(j) + length(loss_inds);
        q = q(survived_inds,:); % take only indices that are left
        
        if(add_new_alleles) % NEW! Add newly born alleles
            num_new_alleles = poissrnd( (block_size /(2*mean_time_allele_polymorphic_at_equilibrium)) * (N_vec(j+1)/N) ); % Proportional to mutation rate times # of chromosomes . Mutation rate is cancelled !
            q = [q' zeros(num_generations, num_new_alleles, 'single')]'; % set frequency for new alleles
            num_simulated_polymorphic_alleles_vec(j+1) = num_simulated_polymorphic_alleles_vec(j+1)+num_new_alleles;
            q((end-num_new_alleles+1:end), j+1) = 1; % set frequency for new alleles
            first_time_vec = [first_time_vec(survived_inds) repmat(j, 1, num_new_alleles)];
            last_time_vec = [last_time_vec(survived_inds) repmat(j, 1, num_new_alleles)];
        end
        
        %%%       weights = weights(survived_inds,:); % take only indices that are left
        %%% q(absorption_inds,j+1) = 1; % /(2*N_vec(j+1));  % start over with new alleles for extint/fixed alleles
        %%% weights(absorption_inds,j+1:end) = N_vec(j+1)/N; % give newly born alleles higher weights
    end % loop on generations
    cur_num_alleles_survived = min(size(q, 1), iters - num_alleles_simulated);
    %     if(size(q,1) > iters  - num_alleles_simulated) % adjust last block
    %         num_simulated_polymorphic_alleles_vec = num_simulated_polymorphic_alleles_vec - size(q,1) ..; % adjust num alleles vec
    %     end
    
    final_q((num_alleles_simulated+1):(num_alleles_simulated+cur_num_alleles_survived),:) = ...
        q(1:cur_num_alleles_survived,:);
    num_alleles_simulated = num_alleles_simulated + cur_num_alleles_survived
    
    j=num_generations; % compute again for last generation
    [U, C] = unique_with_counts(vec2row(q(:,j))); % Compute histogram of counts
    [x_vec{j}, p_vec{j}] = union_with_counts(x_vec{j}, p_vec{j}, [0 U 2*N_vec(j)], [num_losses C num_fixations]);
    % total_het_at_each_generation_vec(j) = 2 .* sum(q(:,j) ./ (2.*N_vec(j)) .* (1-  q(:,j) ./ (2.*N_vec(j)) ) ); % this indicates how much heterozygosity was absorbed at each time
    %     if(~isempty(absorption_inds))
    %         total_het_at_each_generation_vec(j) = total_het_at_each_generation_vec(j) - ...
    %             2 .* sum(q(absorption_inds,j) ./ (2.*N_vec(j)) .* (1-  q(absorption_inds,j) ./ (2.*N_vec(j)) ) ); % this indicates how much heterozygosity was absorbed at each time
    %     end
    ctr_alleles_blocks_simulated = ctr_alleles_blocks_simulated + block_size
end % while num_alleles <= iters

total_het_at_each_generation_vec = total_het_at_each_generation_vec ./ ...
    num_simulated_polymorphic_alleles_vec(1); % normalize by number of iterations
num_absorptions = sum(num_absorptions_by_generation_vec); % count absorptions, losses and fixations
num_fixations = sum(num_fixations_by_generation_vec);
num_losses = sum(num_losses_by_generation_vec);

for j=1:num_generations % sort and normalize distributions
    [x_vec{j}, sort_perm] = sort(x_vec{j});
    p_vec{j} = p_vec{j}(sort_perm) ./ num_simulated_polymorphic_alleles_vec(1); % normalize to incorporate prob. of alleles which started as polymorphic
end
q=final_q; % Copy results

switch init_str % Normalize
    case 'equilibrium'
        total_het_at_each_generation_vec = total_het_at_each_generation_vec .* prob_site_polymorphic_at_equilibrium; % normalize to include fraction of polymorphic sites at the beginning equilibrium
    case 'newly_born'
        total_het_at_each_generation_vec = total_het_at_each_generation_vec .* 2*N*mu; % normalize to include fraction of polymorphic sites born at each generation
        total_het_at_each_generation_vec = cumsum(total_het_at_each_generation_vec); % unite the contribution from each generation
end
prob_site_polymorphic_at_end = prob_site_polymorphic_at_equilibrium * ... esitmate Prob. poly
    num_simulated_polymorphic_alleles_vec(end) / num_simulated_polymorphic_alleles_vec(end);
absorption_time_given_init_freq_vec = absorption_time_given_init_freq_vec ./ count_vec; % normalize to get estimated absorption times


% Compute allele frequency distribution numerically.
% Method: for each generation propagate the markov chain
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
% compute_matrix - flags saying if .. ???
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
function [x_vec, p_vec, total_het_at_each_generation_vec, ...
    num_absorptions, num_fixations, num_absorptions_by_generation_vec, count_vec, ...
    absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, ...
    num_losses, loss_time_given_init_freq_vec, frac_polymorphic_vec, prob_site_polymorphic_at_end, M] = ...
    compute_numeric_expansion_internal(s, mu, two_side_flag, N_vec, init_str, compute_matrix)

num_absorptions = 0; num_fixations = 0; num_losses = 0; % temp
num_absorptions_by_generation_vec = []; count_vec = [];

N = N_vec(1); num_generations = length(N_vec)-1; max_N = max(N_vec);

x_vec = cell(num_generations+1, 1); p_vec = x_vec; het_vec = x_vec;  % structures for keeping track of alleles
total_het_at_each_generation_vec = zeros(num_generations, 1, 'single');
x_vec{1} = (0:2*N) ./ (2*N); % vector of allele frequencies (include monomorphic alleles)
M=[]; 

init_p_vec = exp( allele_freq_spectrum(x_vec{1}, s, N, two_side_flag, 'log') ); % get stationary distrubution at equilibrium (initial condition)
switch init_str
    case 'newly_born'
        init_p_vec(:) = 0; init_p_vec(1) = 1; % Set allele freq. at 0. Contribution will come from 1/2N
    case 'equilibrium'
        init_p_vec(end) = 0; % 1 allele frequency
        init_p_vec = init_p_vec .* 2.* max(10^(-10), mu);  % multiply by theta. (If mu=0: no mutations, then we've got a problem ... set some mu>0)
        init_p_vec(1) = 1 - sum(init_p_vec(2:end-1)); % set first value to complete distribution to sum to one
end % switch init_str

new_x_mean_vec = x_vec{1} .* (1+s); % first get the mean vector of new #offspring 
num_sigmas = 5; % keep only these many st.d. (to speed-up computation)
slack_generations = 0; % burn-in steps
std_vec = sqrt(2*N .* new_x_mean_vec .* (1-new_x_mean_vec)); % binomial standard deviation 
left_range_vec = max(1, round(new_x_mean_vec.*2*N - num_sigmas .* std_vec));
right_range_vec = 1+min(2*N, round(new_x_mean_vec.*2*N + num_sigmas .* std_vec)); % +/- 1
for j=1:slack_generations % perform a few 'burn-in' steps at constant size to reach discrete equilibrium
    run_slack_generation = j
    
    new_p_vec = zeros(2*N+1,1); %
    for k=1:2*N+1 % -1 % loop on current allele (heavy loop)
        %                k_is = k
        cur_range_vec = left_range_vec(k):right_range_vec(k);
        new_p_vec(cur_range_vec) = ...
            new_p_vec(cur_range_vec) + init_p_vec(k)  .* ...
            binopdf(cur_range_vec, 2*N, new_x_mean_vec(k))';
    end
    new_p_vec = vec2column(new_p_vec);
    new_p_vec(1) = new_p_vec(1) + new_p_vec(end); % identify zero and one frequencies
    new_p_vec(end) = 0;
    new_p_vec(2) = new_p_vec(2) + new_p_vec(1) .* 2*N*mu; % correct mutations
    new_p_vec(1) = new_p_vec(1) .* (1-2*N*mu);
    init_p_vec = new_p_vec;
end % slack generations
if(~exist('mu', 'var') || isempty(mu))
    mu =  sum(vec2row(init_p_vec) .* binopdf(0, 2*N, x_vec) ./ (2*N)); % compute mutation rate per site (assuming T=1 target size)
end
p_vec{1} = vec2column(init_p_vec); het_vec{1} = 2.* p_vec{1} .* vec2column(x_vec{1} .* (1-x_vec{1}));
normpdf_vec = normpdf(-6:10^(-5):6); % Prepare in advance Gaussian density to save time. Current resolution: 1 million
for j=1:num_generations % loop on # of generations
    t_gen = cputime;
    p_vec{j} = p_vec{j} ./ sum(p_vec{j}); % normalize to sum to one. Why?
    if(mod(j, ceil(500/N)) == 0)
        run_generation = j
    end
    %            het_vec{j} = normalize_hist(vec2column(x_vec), het_vec{j}); %            Don't normalize !!!!
    x_vec{j+1} = (0:2*N_vec(j+1)) ./ (2*N_vec(j+1)); % new: include also 0 and 1 freqs. to get absolute values
    new_x_mean_vec = x_vec{j} .* (1+s); % first get the mean vector for next generation
    std_vec = sqrt(2*N_vec(j+1) .* new_x_mean_vec .* (1-new_x_mean_vec));
    p_vec{j+1}  = zeros(1, 2*N_vec(j+1)+1); 
    left_range_vec = max(1, round(new_x_mean_vec.*2*N_vec(j+1) - num_sigmas .* std_vec)); % calculate range to compute binomial distribution 
    right_range_vec = 1+min(2*N_vec(j+1), round(new_x_mean_vec.*2*N_vec(j+1) + num_sigmas .* std_vec)); 
    if(compute_matrix)
        M = zeros(2*N_vec(j+1)+1, 2*N_vec(j)+1);
    end
    
    non_negligile_x_inds = find(p_vec{j} > 10^(-9));
    
    skip_running_fraction_of_inds = (2*N_vec(j)+1-length(non_negligile_x_inds)) / (2*N_vec(j)+1) % see how much time we've saved 
    max_range = max(right_range_vec - left_range_vec)
    prob_compute_method = 'approximate'; %'exact'; % s'exact'; % ''; % 'exact'; 
    for k=vec2row(non_negligile_x_inds) % 1:2*N_vec(j)+1 % -1 % loop on current allele frequency. This is the heaviest part
        %        k_is = k
        cur_range_vec = left_range_vec(k):right_range_vec(k);
        
        % Exact binomial computation
        switch prob_compute_method
            case 'exact'
                p_vec{j+1}(cur_range_vec) = ...
                    p_vec{j+1}(cur_range_vec) + p_vec{j}(k)  .* ...
                    binopdf(cur_range_vec-1, 2*N_vec(j+1), new_x_mean_vec(k));
            case {'approx', 'approximate'}
                if(new_x_mean_vec(k) * 2*N_vec(j+1) < 50)  % use poisson approximation  to binomial on left tail
                    p_vec{j+1}(cur_range_vec) = ...
                        p_vec{j+1}(cur_range_vec) + p_vec{j}(k)  .* ...
                        poisspdf(cur_range_vec-1, 2*N_vec(j+1) * new_x_mean_vec(k));  % This part is slowest (binomial). We replace with poisson approximation for larger values!!!
                else % for larger values use poisson or Gaussian approximations
                    if((1-new_x_mean_vec(k)) .* 2*N_vec(j+1) < 50) %  use poisson approximation to binomial on right tail,
                        p_vec{j+1}(cur_range_vec) = ...
                            p_vec{j+1}(cur_range_vec) + p_vec{j}(k)  .* ...
                            poisspdf(2*N_vec(j+1) - (cur_range_vec-1), 2*N_vec(j+1) * (1-new_x_mean_vec(k)));  % This part is slowest. We replaced with poisson !!!
                    else % use Gaussian approximation to binomial in the middle of distribution
                        p_vec{j+1}(cur_range_vec) = ...
                            p_vec{j+1}(cur_range_vec) + p_vec{j}(k)  .* ...
                            normpdf_vec( round(( 6+(cur_range_vec-1 - 2*N_vec(j+1) * new_x_mean_vec(k)) ./ std_vec(k) ) .* 10^5) ) ./ std_vec(k);
                        %                 normpdf(cur_range_vec-1, 2*N_vec(j+1) * new_x_mean_vec(k), ...
                        %                 sqrt(2*N_vec(j+1) * new_x_mean_vec(k) * (1-new_x_mean_vec(k))) );  % This part is slowest. We replace with poisson !!!
                    end
                end % if: what probability propagating approximation to use
        end % switch approximate
        
        if(compute_matrix)
            M(k,:) = binopdf(0:2*N_vec(j+1), 2*N_vec(j+1), new_x_mean_vec(k)); % Why column here? % -1
        end
    end % loop on k
    p_vec{j+1}(1) = p_vec{j+1}(1) + p_vec{j+1}(end); p_vec{j+1}(end) = 0; % combine zero and one frequencies        
    p_vec{j+1}(2) = p_vec{j+1}(2) + p_vec{j+1}(1) .* 2*N_vec(j+1)*mu; % add mutations
    p_vec{j+1}(1) = p_vec{j+1}(1) .* (1-2*N_vec(j+1)*mu);
    if(compute_matrix)
        M(:,1) = M(:,1) + M(:, 2*N_vec(j+1)+1);
        M = M(1:(2*N_vec(j+1)), 1:(2*N_vec(j+1))); % identify ancestral and derived !!! 
        M(:,2) = M(:,2) + M(:,1) .* 2*N_vec(j+1)*mu; % correct for mutation
        M(:,1) = M(:,1) .* (1-2*N_vec(j+1)*mu);
        T_mean = inv(eye(2*N_vec(j+1)-1)-M(2:end,2:end)) * ones(2*N_vec(j+1)-1,1); % caclulate mean absorption time for each frequency
        [V, LAMBDA] = eig(M');  V = V(:,1) ./ sum(V(:,1)); P_POLY = 1-V(1)
        FRAC_LOST = M(2:end,1)'*V(2:end) ./ sum(V(2:end))
    end
    p_vec{j+1} = vec2column(p_vec{j+1});
    p_vec{j+1}(1) = p_vec{j+1}(1) + ... % normalization
        (N_vec(j+1)/N) * (2*N*mu); % add mass at 1/2N frequency. New! add mass to get total frequency above 1
    p_vec{j+1} = p_vec{j+1} ./ sum(p_vec{j+1});
    
    het_vec{j+1} = 2.* p_vec{j+1} .* vec2column(x_vec{j+1} .* (1-x_vec{j+1})); % update heterozygosity     
%    total_het_at_each_generation_vec(j) = sum(het_vec{j}); % 2 .* sum(q(:,j) ./ (2.*N_vec(j)) .* (1-  q(:,j) ./ (2.*N_vec(j)) ) ); % this indicates how much heterozygosity was absorbed at each time
    generation_time = cputime - t_gen
end % loop on generations

p_vec_equlibrium_analytic = p_vec{1};

% Compute summary statistics
total_het_at_each_generation_vec  = sum_cell(het_vec); % total_het_at_each_generation_vec(j+1) = sum(het_vec{j+1});

frac_polymorphic_vec = zeros(num_generations,1);
for j=1:num_generations
    frac_polymorphic_vec(j) = sum(p_vec{j});
end
if(compute_matrix)
    M2 = M; M2([1 end],:) = 0; ones_vec = ones(2*N,1); ones_vec(1) = 0; ones_vec(end) = 0;
    %        ttt=cputime; inv_M2 = inv(M2-eye(2*N+1)); absorption_time_given_init_freq_vec = -sum( inv_M2(:,2:end-1),2); cputime - ttt
    ttt=cputime; absorption_time_given_init_freq_vec = linsolve( M2-eye(2*N),  -ones_vec); inverse_time = cputime - ttt
else
    absorption_time_given_init_freq_vec = zeros(2*N+1,1);
end
prob_fixation = 1; % THIS IS WRONG. HOW TO COMPUTE THIS?
fixation_time_given_init_freq_vec = zeros(size(absorption_time_given_init_freq_vec)); 
loss_time_given_init_freq_vec = zeros(size(absorption_time_given_init_freq_vec)); % Time conditioned on fixation/loss. NEED TO FILL THESE

for j=1:num_generations+1  %change output format (convention)
    p_vec{j} = vec2row(p_vec{j});
    x_vec{j} = round(x_vec{j} .* (2*N_vec(j))); % get number of individual alleles (not fracion)
end
prob_site_polymorphic_at_end = []; % NEED TO FILL !!! 




% Compute everything analytically (how? look at formulas from Shamil for moments)
% NOT WORKING YET !!!
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
% % % % function [total_het_at_each_generation_vec, num_absorptions, num_fixations, ...
% % % %     num_absorptions_by_generation_vec, count_vec, ...
% % % %     absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, num_losses, ...
% % % %     loss_time_given_init_freq_vec, num_simulated_polymorphic_alleles_vec] = ...




% function  [q, weights, x_vec, p_vec, total_het_at_each_generation_vec, ...
%     num_absorptions, num_fixations, num_losses, ...
%     num_absorptions_by_generation_vec, num_fixations_by_generation_vec, num_losses_by_generation_vec, ...
%     absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, loss_time_given_init_freq_vec, ...
%     num_simulated_polymorphic_alleles_vec, count_vec] = ...
function  [x_vec, p_vec, total_het_at_each_generation_vec, num_absorptions ,num_fixations, ...
    num_absorptions_by_generation_vec, count_vec ...
    absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, num_losses, ...
    loss_time_given_init_freq_vec, frac_polymorphic_vec, prob_site_polymorphic_at_end] = ...
    compute_analytic_expansion_internal( ...
    s, mu, two_side_flag, num_generations, N_vec, init_str, num_moments)

use_moments = 1;
use_gegenbauer = 0;
N = N_vec(1); max_N = max(N_vec); 

[total_het_at_each_generation_vec, num_absorptions ,num_fixations, ...
    num_absorptions_by_generation_vec, count_vec ...
    absorption_time_given_init_freq_vec, fixation_time_given_init_freq_vec, num_losses, ...
    loss_time_given_init_freq_vec, frac_polymorphic_vec] = deal([]);
x_vec = cell(num_generations, 1); p_vec = x_vec; 
prob_site_polymorphic_at_equilibrium = (2*N*mu) * 2 * ...
    absorption_time_by_selection(abs(s), 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start


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
    prob_site_polymorphic_at_end = prob_site_polymorphic_at_equilibrium; % initialize
    if(~exist('num_moments', 'var') || isempty(num_moments))
        num_moments = 5; % how many moments to compute
    end
    for j=1:num_generations
        run_j = j
        x_vec{j} = (1:(2*N_vec(j)-1)) ./ (2*N_vec(j));
        [mu_vec_expansion_analytic, mu_vec_equilibrium] = ...
            FisherWright_Compute_SFS_Moments(N_vec(1:j), 0, num_moments); % compute moments with Formulas from Ewens
        % Estimate density using the max-entropy method:
        if(j==1)
            lambda0 = [];
        else % 'warm' start: take initial condition from previous generation
            lambda0 = lambda_max_ent;
        end
        [lambda_max_ent, g_het_max_ent, entr_max_ent] = ...  % Fit density using moments % , lambda0)
            me_dens2(mu_vec_expansion_analytic(2:end) ./ mu_vec_expansion_analytic(1), ...
            x_vec{j}, lambda0, 0); % normalize by 0th moment. Don't plot anything
        
        f_max_ent = zeros(1, 2*N_vec(j)-1);
        for k=1:num_moments
            f_max_ent = f_max_ent + lambda_max_ent(k) .* x_vec{j} .^ (k-1);
        end
        f_max_ent = exp(-f_max_ent) ./ (x_vec{j} .* (1-x_vec{j})); % Compute f. How to normalize?
        f_max_ent = f_max_ent ./ sum(f_max_ent);         
        p_absorb = sum ( f_max_ent .* (x_vec{j}.^(2*N_vec(j)) + (1-x_vec{j}).^(2*N_vec(j))) );
        prob_site_polymorphic_at_end = prob_site_polymorphic_at_end * (1-p_absorb) + ... 
            (1-prob_site_polymorphic_at_end) * 2*N_vec(j)*mu; % XXX UPDATE PROB.        
        p_vec{j} = [0 f_max_ent 0]; % get p. Normalize and add zeros at boundaries 
        x_vec{j} = [0 round(x_vec{j}*2*N) 2*N]; % get x in integers        

    end % loop on number of generations
end % use moments. Function doesn't use mutation rate mu !! 


% Temp: Plot individual trajectories of simulations
%
% Input:
% N - population size
% mu - mutation rate
% s - selection coefficient
% q - individual allele-frequency simulations
% p_vec - probability at each allele frequency
% p_vec_equlibrium_analytic - theoretical solution for population at equilibrium
%
function plot_expansion_internal(N, mu, s, q, x_vec, p_vec, p_vec_equlibrium_analytic)

params_str = ['N=' num2str(N) ', \mu=' num2str(mu) ', s=' num2str(s)];
figure; hold on; plot(q(1:20,:)'); xlabel('Time (generations)'); ylabel('Freq.');
title(['Trajectorues of Fisher-Wright model with selection and mutations. ' params_str]);
figure;  imagesc(q); colorbar;
title(['Trajectorues of Fisher-Wright model with selection and mutations. ' params_str]);
xlabel('Time (generations)'); ylabel('Iteration');

% Need to set x_vec:

figure; hold on;
plot(x_vec(2:(end-1)), p_vec);
plot(x_vec(2:(end-1)), p_vec_equlibrium_analytic, 'r');
xlabel('derived allele freq.'); ylabel('density');
legend('empirical', 'theoretical');
title(['Density of alleles for each allele frequency. ' params_str ...
    ' (' num2str(num_generations) ' generations, ' num2str(iters) ...
    ' iterations, ' num2str(num_effective_iters) ' effective iters)']);

























