% Fit a simple demographic model to population allele frequency spectrum
% We assume that observed alleles are not under selection (neutral)
%
% Input:
% k_vec - number of dervied allele carriers observed in population for each allele
% n_vec - number of individuals profiled in population for each allele
% weights_vec - weight of each allele (default: 1)
% mu - total mutation rate for region (assume this is known for now)
% L_correction_factor - correction factor for empirical heterozygosity
% D_opt - optimization method
%
% Output:
% D - demographic model (expansion, population size etc.)
% max_LL - log-likelihood of data for spectrum
% N_vec_hat - best demographic model (MLE)
%
function [D, max_LL, N_vec_hat, log_like_mat] = fit_demographic_parameters_from_allele_spectrum( ...
    k_vec, n_vec, weights_vec, mu, L_correction_factor, D_opt)

AssignGeneralConstants; AssignRVASConstants;

if(~exist('weights_vec', 'var') || isempty(weights_vec))
    weights_vec = 1;
end
if(length(weights_vec) == length(n_vec))
    weights_vec = weights_vec(n_vec>0); 
end
k_vec = k_vec(n_vec>0); n_vec = n_vec(n_vec>0); % take only polymorphic alleles 

if(~exist('mu', 'var') || isempty(mu))
    mu = mu_per_site; % default mutation rate (per-nucleotide per-generation)
end
if(~exist('L_correction_factor', 'var') || isempty(L_correction_factor))
    L_correction_factor = 1; % ratio between observed and expected num. of polymorphic alleles
end
if(~exist('D_opt', 'var') || isempty(D_opt))
   D_opt = []; D_opt.method = 'forward_simulations';  
end


% We fit expansion model from allele-frequency data, assuming that all alleles are neutral
s = 0; % Assume no selection (synonymous)
alpha=0; % No mixture. Everything is neutral (synonymous)
beta=[];  % effect size on phenotype IRRELEVANT! (we use only genpotypes)

X = [vec2row(k_vec) vec2row(n_vec)]'; % Represent counts in a packed form

D = generate_candidate_models_internal(); % set which candidate models to enumerate 
compute_flag = 'simulation'; % numeric
log_like_mat = zeros(size(D));

filter_by_moments=1; num_moments=1; % First filter by moments !!!
if(filter_by_moments)
    moments_time=cputime;
    het_moment_mat_all_models = zeros(D.num_params, num_moments);
    init_N = zeros(D.num_params, 1);
    for i=1:D.num_params
        if(mod(i, 100) == 0)
            sprintf('Fitting %ld out of %ld demographic model', i, D.num_params)
            cur_moment_time = cputime - moments_time
        end
        D.index = i;
        i_vec = myind2sub(D.num_params_vec, length(D.num_params_vec), i); % create vector of indices
        if(i_vec(2*D.num_stages) == length(D.generations_vec{end})) % get highest number of generations
            j_vec = i_vec;
            for j=length(D.generations_vec{end}):-1:1 % loop on all smaller numbers of generations !
                j_vec(2*D.num_stages) = j;
                cur_ind = mysub2ind(D.num_params_vec, length(D.num_params_vec), j_vec);
                N_vec = demographic_parameters_to_n_vec(D, cur_ind); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation
                if(~filter_N_vec_internal(N_vec)) % get rid of this model
                    max_j = j;
                    break;
                end
            end
            %             if(filter_N_vec_internal(N_vec)) % get rid of this model
            %                 continue;  % PROBLEM! this gets rid of a 'bad' model, and then will get rid of all sub-models, some of which might be 'good'
            %             end
            mu_vec_expansion_analytic = FisherWright_Compute_SFS_Moments(N_vec, 0, num_moments); % compute moments with Formulas from Ewens
            %            het_moment_mat_all_models(i,:) = mu_vec_expansion_analytic(:,end);
            for j=1:max_j % length(D.generations_vec{end}) % loop on all smaller numbers of generations !
                j_vec(2*D.num_stages) = j;
                cur_ind = mysub2ind(D.num_params_vec, length(D.num_params_vec), j_vec);
                het_moment_mat_all_models(cur_ind,:) = mu_vec_expansion_analytic(:,D.generations_vec{3}(j)+D.generations_vec{1}+D.generations_vec{2});
                init_N(cur_ind) = N_vec(1); % take first value of N
            end
        end
    end
    moments_time=cputime - moments_time
end % if filter by moments 

% need to convert sample to population!
%[~, p_vec] = unique_with_counts(k_vec ./ n_vec); % get an empirical distribution
%p_vec = p_vec ./ sum(p_vec); % normalize
% het_moment_mat_data = moment_hist(x_vec, x_vec .* (1-x_vec) .* p_vec, 0, 0, 0) .* L_correction_factor;
het_moment_mat_data = sum ( (k_vec./n_vec) .* (1-k_vec./n_vec) ) * L_correction_factor;  % compute empirical heterozygosity (corrected)
if(~isfield(D_opt, 'moments_epsilon'))
    D_opt.moments_epsilon = 0.1; % allow error in moments
end
region_het_moment_mat_all_models = het_moment_mat_all_models(:,1) .* 4 .* init_N .* mu;
good_inds = find( (abs(region_het_moment_mat_all_models - het_moment_mat_data(:,1)) ./ het_moment_mat_data(:,1))  < D_opt.moments_epsilon );


i_ctr=1; log_like_mat = zeros(D.num_params, 1)-Inf; compute_time=zeros(D.num_params, 1);
test_model_ctr = 0;
if(isempty(good_inds))
    error('Error!! No Demographic Model Fits SFS Moments!!'); 
end
for i=vec2row(good_inds) % 1:D.num_params
    D.index = i; % dset index of current model: should allow multiple indices!!! 
    N_vec = demographic_parameters_to_n_vec(D, D.index); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation
    
    if(filter_N_vec_internal(N_vec))
        continue;
    end
    test_model_ctr = test_model_ctr+1
    %end
    % Compute likelihood. This is trivial one-class likelihood (no mixture ) so should be fast !!!
    rare_cumulative_per_gene = []; % set dummy variables
    target_size_by_class_vec = [mu, mu, mu]; % [neutral, null, missense]
    full_flag = 0; % use summary statistics
    null_w_vec = 1; % NULL_C % assume all alleles are 'null' (but s=0 so actually neutral)
    
    loglike_params = struct('null_w_vec', null_w_vec, 'include_phenotype', 0, ...
            'full_flag', full_flag, 'num_individuals', []);    
    [log_like_mat(i), ~, compute_time(i)] = ... % compute likelihood (here vary only alpha)
        compute_two_class_log_likelihood(s, alpha, beta, target_size_by_class_vec, D, ...
        X, [], [], loglike_params); % null_w_vec, 0, full_flag, []); % don't include phenotype !!
    fprintf('run good ind %ld out of %ld, ', i_ctr, length(good_inds));
    fprintf(' Loglike=%f, cur-time=%f, total-time=%f\n', log_like_mat(i), compute_time(i), sum(compute_time(1:i)));
    
    i_ctr=i_ctr+1;
end

[max_LL, max_J] = max(log_like_mat(good_inds)); max_J = good_inds(max_J); % maximize likelihood. (Take only good inds)
% D = D{max_J};

D.index = max_J;
N_vec_hat = demographic_parameters_to_n_vec(D, D.index);


% % New: allow D_opt to use the fastneutrino method of Song et al.:
% switch D_opt
%     case 'neutrino'
%         R_fastNeutrino = ; % prepare output file
%         savecellfile(R_fastNeutrino, tmp_file);
%         fit_run_str = ['fastNeutrino --maxB 20 --modelFile model1ExponentialEpoch.txt ' ...
%           '--inferredModelOutputFile
%           model1ExponentialEpochInferredParams.txt --maxRandomRestarts 10
%           < sfsNelsonEtAl.txt'];
%
% 'fastNeutrino'; % ...; % generate running command
% end


% Internal function for filtering out obviously unrealistic demographic models.
% Models with N_vec too large or with too many generations are filtered out
% Input:
% N_vec - vector of population sizes
%
% Output:
% filter_flag - 1: filter out this model, 0: keep it
%
function filter_flag = filter_N_vec_internal(N_vec)    % Filter first unreasonable models !!!
filter_flag=0;
if(max(N_vec) > 10^9) % over 1 billion individuals 
    filter_flag=1;
end
if(max(N_vec) > 10^4 * N_vec(1)) % don't allow too big expansion
    filter_flag=1;
end
if(length(N_vec) > 4000) % we model maximum of 4000 generations
    filter_flag=1;
end


% Internal function for generating possible demographic models 
% 
% Output: 
% D - structure with candidate models 
function D = generate_candidate_models_internal()

% Set different demographic models:
D.init_pop_size_vec{1} = 500; % 1. ancestral population
D.init_pop_size_vec{2} = round(logspace(1, 3, 7)); % 2. bottleneck size 1, 4, 10
D.init_pop_size_vec{3} = -1; % 3. population after bottleneck
D.init_pop_size_vec{4} = -1; % 3. population in second expansion phase 


D.generations_vec{1} = 500; % 4. burnout time
D.generations_vec{2} = 10;  % 5. number of generations since bottleneck
D.generations_vec{3} = round(logspace(1, 4, 10));  % 6. number of generations for expansion
D.generations_vec{4} = 100;

D.expan_rate_vec{1} = 1; % fixed ancestral population size
D.expan_rate_vec{2} = 1; % fixed buttleneck size
D.expan_rate_vec{3} = logspace(0.001, 0.033333, 5); % 5. expansion rate per generation
D.expan_rate_vec{4} = 1; % 5. expansion rate per generation in second phase 

D.num_params_vec = [length_cell(D.init_pop_size_vec) length_cell(D.generations_vec) length_cell(D.expan_rate_vec)];
D.num_stages = length(D.init_pop_size_vec);

D.num_params = prod(D.num_params_vec); % how many parameters to enumerate
D.use_allele_counts = 0; % new! don't use counts in likelihood computations !!!

