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
D.generations_vec{3} = round(logspace(1, 4, 20));  % 6. number of generations for expansion
%D.num_params_vec(2) = length(D.generations_vec{3});

D.expan_rate_vec{1} = 1; % fixed ancestral population size
D.expan_rate_vec{2} = 1; % fixed buttleneck size
D.expan_rate_vec{3} = logspace(0.001, 0.1, 20); % 5. expansion rate per generation
%D.num_params_vec(3) = length(D.expan_rate_vec{3});

D.num_params_vec = [length_cell(D.init_pop_size_vec) length_cell(D.generations_vec) length_cell(D.expan_rate_vec)];
D.num_stages = length(D.init_pop_size_vec);

%D.num_bottleneck_size = length(D.bottleneck_size);
%D.num_bottleneck_size = length(D.bottleneck_size);
%D.num_init_pop_size = length(D.bottleneck_size);
%D.num_expansion_rate = length(D.expansion_rate);


D.num_params = prod(D.num_params_vec); % how many parameters to enumerate
compute_flag = 'numeric';

log_like_mat = zeros(size(D));

filter_by_moments=1; num_moments=1; % First filter by moments !!!
if(filter_by_moments)
    moments_time=cputime;
    het_moment_mat_all_models = zeros(D.num_params, num_moments);
    for i=1:D.num_params
        if(mod(i, 100) == 0)
            sprintf('Fitting %ld out of %ld demographic model', i, D.num_params)
            cur_moment_time = cputime - moments_time
        end
        D.index = i;
        i_vec = myind2sub(D.num_params_vec, length(D.num_params_vec), i); % create vector of indices
        if(i_vec(2*D.num_stages) == length(D.generations_vec{end})) % get highest number of generations
            N_vec = demographic_parameters_to_n_vec(D, D.index); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation
            [mu_vec_expansion_analytic] = FisherWright_Compute_SFS_Moments(N_vec, 0, num_moments); % compute moments with Formulas from Ewens
            het_moment_mat_all_models(i,:) = mu_vec_expansion_analytic(:,end);
            
            for j=1:length(D.generations_vec{end}) % loop on all smaller numbers of generations !
                cur_ind=j; % TEMP WRONG !!
                het_moment_mat_all_models(cur_ind,:) = mu_vec_expansion_analytic(:,D.generations_vec{3}(j));
            end
            
        end
    end
    moments_time=cputime - moments_time
end


% need to convert sample to population!
het_moment_mat_data = moment_hist(x_vec, x_vec .* (1-x_vec) .* p_vec, 0, 0, 0);



epsilon = 0.05; % allow error in moments
good_inds = find( abs(het_moment_mat_all_models(:,1) - het_moment_mat_data(:,1))  < epsilon );

for i=good_inds % 1:D.num_params
    
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
    
    [demographic_x_vec, demographic_f_vec] = compute_allele_freq_spectrum_from_demographic_model(D, 0, compute_flag); % Try a grid of different values
    
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




