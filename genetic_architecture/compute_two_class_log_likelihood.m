% Compute log-likelihood for data using the two-class model
% Formula is taken from eq. (XXX) from document on RVAS:
%
% Input:
% s_null_vec - vector of possible selection coefficients for the null alleles. Should be NEGATIVE for deleterious alleles (so fitness is 1+s)
% alpha_vec - vector of fraction of null alleles at birth
% beta_vec - vector of effect size of null alleles
% target_size_by_class_vec - NEW! target size for each allele type (used in poisson model).
%                            This is like the above but counts expected NUMBER of different alleles,
%                            not their FREQUENCIES !  %%% NEW! alleles class type vector. This states for each allele if it is neutral, harmfull, or in a mixture (missense)
% D - demographic model (could be N:  effective population size for equilibrium model)  NEW! (currently only constant pop. size supported)
% X - genotype data matrix - for some models we actually need it. For others, just the allele frequencies.
% y - phenotype data vector - (optional) may be used if likelihood includes phenotype
% trait_struct - structure with: type - disease or quantitative
%                                prevalence - frequency of disease in the population for disease trait
% null_w_vec - (optional) this is the assignment of which alleles are null and which not.
%                   This simplifies the likelihood computation a lot when it is known
%                   (no need to average over possible assignments of w). Convention:
%                      1: null
%                      0: neutral
%                     -1: missense (unknown, drawn from a mixture)
% include_phenotype - flag saying which part of the likelihood we compute:
%                      1: all (genotype+phenotype) (default)
%                      0: only genotype part of LogLikelihood
%                     -1: only phenotype part of LogLikelihood
%                       : NEW !!! Allow flag for only genotype, but without
%                       counting (only SFS!!)
% full_flag - input format:
%                      1: X is genotype matrix (default)
%                      0: X contains sufficient statistics (sums of rows and columns)
% num_individuals - input number of individuals in sample (when not given in input vectors)
% print_flag - print likelihood (optional)
%
% Output:
% log_like_mat - Matrix (3-d) of log-likelihood of data for each parameter choice (s, alpha and beta)
% P_poly - Structure with computation information. Probability of polymorphic alleles and more information
% compute_time - total time operation took
%
function [log_like_mat, P_poly, compute_time] = ...
    compute_two_class_log_likelihood(s_null_vec, alpha_vec, beta_vec, ...
    target_size_by_class_vec, D, ...
    X, y, trait_struct, params)  % null_w_vec, include_phenotype, full_flag, num_individuals, print_flag)

%%%%%%%%%%%%%%%%%%%%%%%% Set Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttt=cputime;
use_allele_freq_flag=2; % should be 2. (temp for debugging - allow computing only partial likelihoods)
AssignGeneralConstants; AssignRVASConstants;

params = internal_set_parmas(params);

if(~isfield(trait_struct, 'type'))
    trait_struct.type = 'quantitative';
end
if(~isfield(trait_struct, 'prevalence'))
    trait_struct.prevalence = [];
end
if(~exist('target_size_by_class_vec', 'var') || isempty(target_size_by_class_vec))
    poisson_model_flag = 0; expand_format_flag = 'individual'; % include individual-level information (how many alleles for each individual)
else
    poisson_model_flag = 1; expand_format_flag = 'summary'; % use poisson model for each class
end
if(isstruct(D)) % Set compute method
    compute_flag = 'simulation'; % simulate forward alleles
else % equilibrium model
    N=D; compute_flag = 'analytic';
    D = []; D.N_vec = N;
end
if(~isfield(D, 'use_allele_counts'))
    D.use_allele_counts = 1;
end
if(~isfield(D, 'mu'))
    D.mu = mu_per_site; % TEMP! Set a default value for mu !
end
if(~params.include_phenotype) % here beta doesn't influence result
    beta_vec = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%% End Set Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_s = length(s_null_vec);
num_alpha = length(alpha_vec);
num_beta = length(beta_vec);
num_classes = 3; % 1 - synonymous (neutral), 2 - null, 3 - missense ?

if(params.full_flag)
    [params.num_individuals, L] = size(X); % set number of individuals and number of SNPs
else
    if(~isempty(y)) % when phenotype is given
        params.num_individuals = length(y);
    else % here we have no way of knowing the # of individuals !!!
    end
end
if(params.full_flag)
    x_inds = cell(params.num_individuals,1);
    for i=1:params.num_individuals % record alternative alleles
        x_inds{i} = find(X(i,:));
    end
    num_carriers_vec = sum(X); % num. of individuals carrying each rare allele
    num_alleles_vec = sum(X,2); % num. of rare alleles in each individual
    num_individuals_vec = repmat(params.num_individuals, L, 1); % just repeat (assume for all alleles the same # of individuals is profiled)
else % here we only have summary statistics (sum of rows and columns)
    [num_alleles_vec, ~, num_carriers_vec, num_individuals_vec, L] = ...
        expand_two_class_summary_statistics(X, params.num_individuals, expand_format_flag); % get num. of individuals for each rare allele, and num. of rare allele for each individual
end
if(isscalar(params.null_w_vec)) % all alleles with same class
    params.null_w_vec = repmat(params.null_w_vec, L, 1);
end

n =  max(num_individuals_vec);
[unique_num_alleles_vec, ~, J_num_alleles] = unique(num_alleles_vec); 
[unique_num_carriers, I_num_carriers, h_num_carriers] = get_duplicates(num_carriers_vec);
if(min(unique_num_carriers) > 0) % add monomorphic alleles
    params.null_w_vec(end+1) = params.null_w_vec(1);
    unique_num_carriers(end+1) = 0;
    h_num_carriers(end+1) = 0;
    I_num_carriers{end+1} = length(params.null_w_vec);
end
if(max(unique_num_carriers) < n) % add monomorphic alleles
    params.null_w_vec(end+1) = params.null_w_vec(1);
    unique_num_carriers(end+1) = n;
    h_num_carriers(end+1) = 0;
    I_num_carriers{end+1} = length(params.null_w_vec);
end

if(poisson_model_flag) % Compute counts
    num_polymorphic_alleles_observed = zeros(3,1);
    for i=1:num_classes % loop on 3 classes
        num_polymorphic_alleles_observed(i) = sum(X(params.null_w_vec+2 == i,:)>0); % polymorphic in SAMPLE for each classs
    end
    P_poly = [];
    [P_poly.population, P_poly.sample] = deal(zeros(num_s, num_classes));
    [P_poly.missense_population, P_poly.missense_sample] = deal(zeros(num_s, num_alpha));
end
log_like_mat = zeros(num_s, num_alpha, num_beta);

% Determine distribution. For equilibrium use analytic solution. Otherwise, use simulation!!!
[x_vec, log_x_vec, log_one_minus_x_vec, allele_freq_hist] = deal(cell(num_classes, 1)); [sum_allele_freq_hist, L_correction_factor ] = deal(zeros(num_classes, 1));
switch compute_flag
    case 'analytic'
        x_vec{NEUTRAL_C} = (1:2*N-1) ./ (2*N); % vector of allele frequencies
        allele_freq_hist{NEUTRAL_C} = exp(allele_freq_spectrum(x_vec{NEUTRAL_C}, 0, N, 0, 'log')); % allele freq. distribution for neutral alleles. NOT Normalized!
    case 'simulation'
        if(~isfield(D, 'N_vec'))
            D.N_vec = demographic_parameters_to_n_vec(D, D.index); 
        end
        N=D.N_vec(1);
        if(isfield(D, 'SFS')) % here specrum already computed
            allele_freq_hist{NEUTRAL_C} = D.SFS.p_vec(1,:); % take neutral
            x_vec{NEUTRAL_C} = D.SFS.x_vec; % copy x vec
            L_correction_factor(NEUTRAL_C) = D.SFS.L; % take target size
        else
            [x_vec{NEUTRAL_C}, allele_freq_hist{NEUTRAL_C}, L_correction_factor(NEUTRAL_C), demographic_compute_time] = ...
                compute_allele_freq_spectrum_from_demographic_model(D, 0, compute_flag); % Try a grid of different values
            sprintf('Compute neutral spectrum time=%f', demographic_compute_time)
        end
        x_vec{NEUTRAL_C} = x_vec{NEUTRAL_C} ./ (2*D.N_vec(end-1)); % normalize: from counts to allele freq.
        
        % NEW: can use moments !!! (so analytic but for non-constant population size)
end
sum_allele_freq_hist(NEUTRAL_C) = sum(allele_freq_hist{NEUTRAL_C});
[sample_x_vec{NEUTRAL_C}, sample_p_vec{NEUTRAL_C}] = population_to_sample_allele_freq_distribution( ...
    x_vec{NEUTRAL_C}, allele_freq_hist{NEUTRAL_C}, max(num_individuals_vec), round(unique_num_carriers)', 1);

prob_null_given_x = zeros(L,1); % conditional probability of allele being null when we know the frequency
T_0 = sum(allele_freq_hist{NEUTRAL_C}) * (x_vec{NEUTRAL_C}(2)-x_vec{NEUTRAL_C}(1)); %T_0 = absorption_time_by_selection(0, 1, N, 1/(2*N), 1-1/(2*N), 0); % use analytic approximation (not histogtam). Turns out to matter a lot! %%% T_0 = integral_hist(x_vec{NEUTRAL_C}, allele_freq_hist{NEUTRAL_C}); %%% T_0 = absorption_time_by_selection(0, 1, N, 1/(2*N), 1-1/(2*N), 0); % 'freq');
log_x_vec{NEUTRAL_C} = log(x_vec{NEUTRAL_C}); log_one_minus_x_vec{NEUTRAL_C} = log(1-x_vec{NEUTRAL_C});
%tmp_likelihood_one_allele = repmat(BIG_NUM, max(num_individuals_vec), num_classes);

if(params.optimize_alpha) % Here find best alpha value
    if(params.include_phenotype) % prepare binomial tables
        binom_vec = zeros(length(unique_num_alleles_vec), max(unique_num_alleles_vec)+1); % save binomial coefficients in table
    end
end
if(poisson_model_flag)
    lambda_s = zeros(num_s, 1);
    lambda_missense = zeros(num_s, num_alpha);
    % tmp_z0_vec = exp( num_individuals_vec(1) .* (log_one_minus_x_vec{NEUTRAL_C})); %  - log(underflow_correction_q(1))  ); % what's this? multinomial/binomial coefficient?
    P_poly.neutral_population = 4 * N * D.mu * T_0; % prob. allele polymorphic in population
    P_poly.neutral_sample = P_poly.neutral_population .* ...
        (1 - sum( allele_freq_hist{NEUTRAL_C} .* (1-x_vec{NEUTRAL_C}) .^ num_individuals_vec(1) ) / sum(allele_freq_hist{NEUTRAL_C}));
    
    %    tmp_z_vec = cell(num_classes, 1);
    %    tmp_z_vec{NEUTRAL_C} = sparse(length(unique_num_carriers), length(x_vec{NEUTRAL_C})); % create sparse matrix
    %    pos_tmp_z_inds = cell(length(unique_num_carriers), num_classes);
    %%    t2 = cputime; % tmp_ind_vec = []; tmp_val_vec = [];
    %    for j=1:length(unique_num_carriers) % loop on unique
    %        tmp_vec = exp( log_binom(num_individuals_vec(I(j)), num_carriers_vec(I(j))) + ...
    %            num_carriers_vec(I(j)) .* log_x_vec{NEUTRAL_C} + ...
    %            (num_individuals_vec(I(j))-num_carriers_vec(I(j))) .* log_one_minus_x_vec{NEUTRAL_C} ); % assumes num individuals is the same
    %        pos_tmp_z_inds{j, NEUTRAL_C} = find(tmp_vec > 10^(-8));
    %        tmp_z_vec{NEUTRAL_C}(j, pos_tmp_z_inds{j, NEUTRAL_C}) = tmp_vec(pos_tmp_z_inds{j, NEUTRAL_C});
    % end
    %    cputime - t2
end

T_s = zeros(num_s,1); % time an allele spends at polymorphic state
for i_s = 1:num_s % loop on parameters
    switch compute_flag     % Get allele frequency based on spectrum
        case 'analytic'
            x_vec{NULL_C} = x_vec{NEUTRAL_C};
            allele_freq_hist{NULL_C} = exp(allele_freq_spectrum(x_vec{NULL_C}, s_null_vec(i_s), N, 0, 'log')); % allele freq. distribution for null alleles.  NOT Normalized!
        case {'simulations', 'simulation'}
            if(s_null_vec(i_s) == 0) % neutral. already computed, save time
                x_vec{NULL_C} = x_vec{NEUTRAL_C}; allele_freq_hist{NULL_C} = allele_freq_hist{NEUTRAL_C};
                log_x_vec{NULL_C} = log_x_vec{NEUTRAL_C}; log_one_minus_x_vec{NULL_C} = log_one_minus_x_vec{NEUTRAL_C};
            else % here s not 0 (deleterious alleles). Compute spectrume again !!
                if(isfield(D, 'SFS')) % here specrum already computed
                    [~, j_s] = min(abs(s_null_vec(i_s)) - abs(D.s_grid));
                    allele_freq_hist{NULL_C} = D.SFS.p_vec(j_s,:); % take neutral
                    x_vec{NULL_C} = D.SFS.x_vec; % copy x vec
                    L_correction_factor(NULL_C) = D.SFS.L; % take target size
                else % here compute SFS using simulations
                    [x_vec{NULL_C}, allele_freq_hist{NULL_C}, L_correction_factor(NULL_C), demographic_compute_time] = ...
                        compute_allele_freq_spectrum_from_demographic_model(D, s_null_vec(i_s), compute_flag); % Try a grid of different values
                    sprintf('Compute null spectrum time=%f', demographic_compute_time)
                end
                x_vec{NULL_C} = x_vec{NULL_C} ./ (2*D.N_vec(end-1)); % normalize: from counts to allele freq.
                log_x_vec{NULL_C} = log(x_vec{NULL_C}); log_one_minus_x_vec{NULL_C} = log(1-x_vec{NULL_C});
            end % if s==0
    end
    [sample_x_vec{NULL_C}, sample_p_vec{NULL_C}] = population_to_sample_allele_freq_distribution( ...
        x_vec{NULL_C}, allele_freq_hist{NULL_C}, max(num_individuals_vec), round(unique_num_carriers)', 1);
    sum_allele_freq_hist(NULL_C) = sum(allele_freq_hist{NULL_C});
    T_s(i_s) = sum(allele_freq_hist{NULL_C}) * (x_vec{NULL_C}(2)-x_vec{NULL_C}(1));     %%    T_s(i_s) = absorption_time_by_selection(-s_null_vec(i_s), 1, N, 1/(2*N), 1-1/(2*N), 0); % use analytic approximation (not histogtam). Turns out to matter a lot!     %    T_s = integral_hist(x_vec{NULL_C}, allele_freq_hist{NULL_C}); %%% T_s = absorption_time_by_selection(-s_null_vec(i_s), 1, N, 1/(2*N), 1-1/(2*N), 0); % 'freq'); % sure we need to use 'freq' here ???
    
    if(poisson_model_flag) % compute poisson part of likelihood
        % we must assume here we know the number of individuals profiled at each region
        %        tmp_z0_vec = exp( num_individuals_vec(1) .* (log_one_minus_x_vec{NEUTRAL_C})); %  - log(underflow_correction_q(1))  ); % what's this? multinomial/binomial coefficient?
        P_poly.population(i_s, NULL_C) = 4 * N * D.mu * T_s(i_s); % prob. allele polymorphic in population. Is this stuff just for equilibrium?
        P_poly.sample(i_s, NULL_C) =  P_poly.population(i_s, NULL_C) .* ...
            (1 - sum( allele_freq_hist{NULL_C} .* (1-x_vec{NULL_C}) .^ num_individuals_vec(1) ) / sum(allele_freq_hist{NULL_C}));
        %            (1 - integral_hist(x_vec{NULL_C},  allele_freq_hist{NULL_C} .* tmp_z0_vec ) ./ integral_hist(x_vec{NULL_C}, allele_freq_hist{NULL_C}) );
        P_poly.missense_population(i_s,:) = alpha_vec .* P_poly.population(i_s, NULL_C) + ...
            (1-alpha_vec) .* P_poly.neutral_population;
        P_poly.missense_sample(i_s,:) = alpha_vec .* P_poly.sample(i_s, NULL_C) + ...
            (1-alpha_vec) .* P_poly.neutral_sample;
        
        % set lambdas for poisson model
        lambda_s(i_s) = 4 * N * D.mu * target_size_by_class_vec(NULL_C) * T_s(i_s);
        lambda_missense(i_s,:) = 4 * N * D.mu * target_size_by_class_vec(MISSENSE_C) .* ( alpha_vec .* T_s(i_s) + (1-alpha_vec) .* T_0 );
    end % if poisson
    
    for i_alpha = 1:num_alpha
        ttt_one = cputime;
        %        p_null = alpha_vec(i_alpha) * T_s(i_s) / (alpha_vec(i_alpha) * T_s(i_s) + (1-alpha_vec(i_alpha)) * T_0); % compute probability that a given observed polymorphic locus is null
        %        allele_freq_hist{MISSENSE_C} = alpha_vec(i_alpha) .* allele_freq_hist{NULL_C} + (1-alpha_vec(i_alpha)) .* allele_freq_hist{NEUTRAL_C}; % mixture of neutral and null freq. distributions.  NOT Normalized! the weight of neutral and null may be very different!
        [x_vec{MISSENSE_C}, allele_freq_hist{MISSENSE_C}] = sum_hist(x_vec{NULL_C}, alpha_vec(i_alpha) .* allele_freq_hist{NULL_C}, ...
            x_vec{NEUTRAL_C}, (1-alpha_vec(i_alpha)) .* allele_freq_hist{NEUTRAL_C}, 1, 2); % add histograms (with different x values. Normalize to keep histogram sums)
        [sample_x_vec{MISSENSE_C}, sample_p_vec{MISSENSE_C}] = population_to_sample_allele_freq_distribution( ...
            x_vec{MISSENSE_C}, allele_freq_hist{MISSENSE_C}, max(num_individuals_vec), unique_num_carriers', 1);
        sum_allele_freq_hist(MISSENSE_C) = sum(allele_freq_hist{MISSENSE_C});
        if(params.include_phenotype && params.optimize_alpha) % prepare binomial tables
            for i=1:length(unique_num_alleles_vec)
                binom_vec(i,1:unique_num_alleles_vec(i)+1) = ...
                    binopdf(0:unique_num_alleles_vec(i), unique_num_alleles_vec(i), alpha_vec(i_alpha));
            end
        end % if include phenotype
        
        if(poisson_model_flag && D.use_allele_counts) % compute poisson part of likelihood  (observation of polymorphic alleles)
            log_like_mat(i_s, i_alpha, :) = internal_log_like_count_polymorphic_alleles(target_size_by_class_vec, num_polymorphic_alleles_observed, ...
                P_poly, use_allele_freq_flag, i_s, i_alpha);
        else
            log_like_mat(i_s, i_alpha, :) = 0;
        end % if poisson model flag
        for i_beta = 1:num_beta % loop on effect size
            %            tmp_likelihood_one_allele(:) = BIG_NUM; % repmat(BIG_NUM, max(num_individuals_vec), 3);
            %            time_before_looping_on_alleles = cputime-ttt
            P_poly.LL_vec =  h_num_carriers .* log(sample_p_vec{params.null_w_vec(1)+2}); % TEMP !!!! for DEBUG !!!
            P_poly.sample_p_vec = sample_p_vec{params.null_w_vec(1)+2}; % TEMP !!!! for DEBUG !!!
            P_poly.counts_vec = h_num_carriers;  % TEMP: record counts of data
            for j=1:length(unique_num_carriers) % here we go over alleles by frequency. We assume all of the same class (w_null)
                if((unique_num_carriers(j)==0) || (unique_num_carriers(j)==n))
                    continue;
                end
                if(params.include_genotype) %  || optimize_alpha)
                    log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                        params.null_w_vec(I_num_carriers{j}(1)) * h_num_carriers(j) * log(sample_p_vec{params.null_w_vec(I_num_carriers{j}(1))+2}(j)); % why + 2? 1 -> NULL_C=3
                end
            end
            %            for j=1:L    % Compute genotype part. Loop on loci. Compute for each # of carriers only once
            %                if( (num_carriers_vec(j) == 0) || (null_w_vec(j) == 0) )% loop only on polymorphic non-synonmous alleles
            %                    continue;
            %                end
            %
            %                num_carriers_ind = find(num_carriers_vec(j) == unique_num_carriers);
            %                if(include_genotype) %  || optimize_alpha)
            %                    %                    if(tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2 ) == BIG_NUM) % this means we've already computed for this frequency
            %                    %                        if(include_genotype)
            %                    %                            tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2) = ...
            %                    %                                log( max(10^(-100), sum(allele_freq_hist{null_w_vec(j)+2}(pos_tmp_z_inds{num_carriers_ind, null_w_vec(j)+2}) .* ...
            %                    %                                tmp_z_vec{null_w_vec(j)+2}(num_carriers_ind, pos_tmp_z_inds{num_carriers_ind, null_w_vec(j)+2})) / sum_allele_freq_hist(null_w_vec(j)+2)) );
            %                    %                        end
            %                    %                    end % already computed likelihood
            %                    log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
            %                        log(sample_p_vec{null_w_vec(j)+2}(num_carriers_ind));
            %                    %                tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2);
            %                end % include genotypes or optimize alpha
            %                if(mod(j, 1000) == 0)
            %                    run_locus = j
            %                end
            %                if(include_phenotype && optimize_alpha)  % compute conditional probability of allele being null given frequency x.
            %                    prob_null_given_x(j) = p_null * integral_hist(x_vec{NULL_C}, allele_freq_hist{NULL_C} .* tmp_z_vec{null_w_vec(j)+2});
            %                    prob_null_given_x(j) = prob_null_given_x(j) / (prob_null_given_x(j) + ...
            %                        (1-p_null) * integral_hist(x_vec{NEUTRAL_C}, allele_freq_hist{NEUTRAL_C} .* tmp_z_vec{null_w_vec(j)+2}));
            %                end
            %    end % loop on loci
            %            ttt_loop_on_loci = cputime - ttt
            if(params.include_genotype) % here add genotype part of log-likelihood
                %                log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + 0 * sum(log_like_correction); % add correction. Just a constant (who cares)
                if(poisson_model_flag) % here include all alleles classes
                    log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) - 0; %%%%... % add normalization constant
                else % here assume all alleles are missnese
                    log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) - ... % add normalization constant
                        L * log( alpha_vec(i_alpha) * T_s(i_s) + (1-alpha_vec(i_alpha)) * T_0 );
                end
            end % if include genotypes
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Start if include phenotypes (skipped) %%%%%%%%%%%%
            if(params.include_phenotype) % compute phenotypes contribution to likelihood
                log_like_mat(i_s,i_alpha,i_beta) = internal_log_like_pheotype(log_like_mat(i_s,i_alpha,i_beta), ...
                    X, x_inds, y, trait_struct, alpha_vec(i_alpha), beta_vec(i_beta), J_num_alleles, params.num_individuals, num_alleles_vec, optimize_alpha, prob_null_given_x);
            end % if to include phenotypes
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% End if include phenotypes %%%%%%%%%%%%
            if(params.print_flag)
                if( (mod(i_alpha, 50) == 0) || mod(i_s, 50) == 0)
                    run_index_s_alpha_beta = [i_s i_alpha i_beta]
                    time_one_likelihood = cputime-ttt_one
                end
            end
        end % loop on effect size beta
    end % loop on mixture coefficient alpha
end % loop on selection coefficients

compute_time = cputime-ttt;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End Main Computations. Start Internal Functions %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Internal function: return the log-likelihood of data (new one used)
%
% Input:
% target_size_alleles -
% num_polymorphic_alleles_observed -
% P_poly - structure with probability of alleles from different class being polymorphic
% use_allele_freq_flag -
% i_s - index of s, selection coefficient
% i_alpha - index of alpha, frac of missense alleles which are null
%
% Output:
% log_like_vec - vector of log-likelihood for each parameter
%
function log_like_vec = internal_log_like_count_polymorphic_alleles( ...
    target_size_alleles, num_polymorphic_alleles_observed, P_poly, ...
    use_allele_freq_flag, i_s, i_alpha)

AssignRVASConstants;
if( ismember(use_allele_freq_flag, [0, 2]))
    log_like_vec = ...
        (target_size_alleles(NULL_C) - num_polymorphic_alleles_observed(NULL_C)) .* ...
        log(1-P_poly.sample(i_s, NULL_C)) + ...
        (target_size_alleles(NEUTRAL_C) - num_polymorphic_alleles_observed(NEUTRAL_C)) .* ... % Compute log-likelihood for alleles not present
        log(1-P_poly.neutral_sample) + ...
        (target_size_alleles(MISSENSE_C) - num_polymorphic_alleles_observed(MISSENSE_C)) .* ...
        log(1-P_poly.missense_sample(i_s, i_alpha));
    
    log_like_vec = log_like_vec + ... % add log-likelihood for observed alleles
        num_polymorphic_alleles_observed(NULL_C) .* ...
        log(P_poly.sample(i_s, NULL_C)) + ...
        num_polymorphic_alleles_observed(NEUTRAL_C) .* ...
        log(P_poly.neutral_sample) + ...
        num_polymorphic_alleles_observed(MISSENSE_C) .* ...
        log(P_poly.missense_sample(i_s, i_alpha));  % all probs here should be in POPULATION!!! (if we include all parts of likelihood!)
    if(use_allele_freq_flag == 0)
        return;
    end
else
    log_like_vec = 0;
end

if( ismember(use_allele_freq_flag, [1, 2])) % add log of binomial coefficient (?)
    log_like_vec = log_like_vec + ...
        num_polymorphic_alleles_observed(NULL_C) .* ...
        log(P_poly.population(i_s, NULL_C)) + ...
        num_polymorphic_alleles_observed(NEUTRAL_C) .* ...
        log(P_poly.neutral_population) + ...
        num_polymorphic_alleles_observed(MISSENSE_C) .* ...
        log(P_poly.missense_population(i_s, i_alpha));            % all probs here should be in POPULATION!!! (if we include all parts of likelihood!)
    
    log_like_vec = log_like_vec - ...
        num_polymorphic_alleles_observed(NULL_C) .* ...
        log(P_poly.sample(i_s, NULL_C)) - ...
        num_polymorphic_alleles_observed(NEUTRAL_C) .* ...
        log(P_poly.neutral_sample) - ...
        num_polymorphic_alleles_observed(MISSENSE_C) .* ...
        log(P_poly.missense_sample(i_s, i_alpha));            % all probs here should be in POPULATION!!! (if we include all parts of likelihood!)
    %                                 log_like_mat(i_s,i_alpha,:) = 0;
end % if use allele freq flag



% Internal function: return phenotype part of log-likelihood
%
% Input:
% log_like_input, ...
% X -
% x_inds -
% y - vector of phenotype
% trait_struct -
% alpha -
% beta -
% J_num_alleles -
% params.num_individuals -
% num_alleles_vec -
% optimize_alpha -
% prob_null_given_x
%
% Output:
% log_like - log likelihood
%
function log_like = internal_log_like_pheotype(log_like_input, ...
    X, x_inds, y, trait_struct, alpha, beta, J_num_alleles, ...
    num_individuals, num_alleles_vec, optimize_alpha, prob_null_given_x)

sigma_e = 1; % environmental noise level

log_like = log_like_input;
switch trait_struct.type
    case {'binary', 'disease'}
        mean_x = mean(X(:));
        mean_f = mean_x * alpha;
        var_exp = beta.^2 * mean_f * (1-mean_f);
        sigma_e = 1 - var_exp;
end
if(optimize_alpha)       % Perform exponential search on all values of w (which alleles are null and which not). Is this the right formula?
    cond_y_tab = zeros(num_individuals, max(num_alleles_vec)+1); % table. entry (i,j) is Prob. (y(i) | beta*(j-1) nulls)
    for i=1:num_individuals
        %                        for j=1:num_alleles_vec(i)+1
        cond_y_tab(i,1:num_alleles_vec(i)+1) = internal_phenotype_fun( ...
            y(i), beta.*(0:num_alleles_vec(i)), sigma_e, trait_struct);
        %                        end
    end
    if(params.full_flag) % exponential search
        W_mat = my_dec2base( 0:2^max(num_alleles_vec)-1, 2, max(num_alleles_vec)); % The heavy part: exponentially many genotypes
        W_sum_cell = cell(max(num_alleles_vec)+1,1); % what is this?
        for i=0:max(num_alleles_vec)
            W_sum_cell{i+1} = sum(W_mat(1:2^i,1:i),2);
        end
        
        for i=1:num_individuals          % Compute phenotype part. Loop on individuals. Heaviest loop
            %                 i_is = i
            %                 num_alleles_is =num_alleles_vec(i)
            k = num_alleles_vec(i);
            prob_null_mat = repmat(prob_null_given_x(x_inds{i}), 1, 2^k)';
            prob_null_mat = prod( prob_null_mat .^ W_mat(1:2^k,1:k) .* (1-prob_null_mat) .^ (1-W_mat(1:2^k,1:k)), 2);
            log_like = log_like + ...
                log(sum(prob_null_mat .* cond_y_tab(i, W_sum_cell{k+1}+1)')); % why plus one???
        end % loop on individuals
    else % use binomial distribution (linear search instead of exponential) % here X is summary statistics (row and column sum)
        for i=1:num_individuals % loop on individuals (heavy part)
            log_like = log_like + ...
                log( sum( binopdf(0:num_alleles_vec(i), num_alleles_vec(i), alpha) .* ...
                cond_y_tab(i,1:num_alleles_vec(i)+1)  ) ); % sum over binomial probabilities
            log_like = log_like + ...
                log( sum( binom_vec(J_num_alleles(i), 1:num_alleles_vec(i)+1) .* ...
                cond_y_tab(i,1:num_alleles_vec(i)+1)  ) ); % sum over binomial probabilities
        end
    end
else % here assume that we know which alleles are null and which are not
    if(params.full_flag) % here X is a genotype matrix
        weight_vec = X * vec2column(null_w_vec); % vector saying how many null alleles are in each individual
    else % here X is just summary statistics (row and column sums)
        switch expand_format_flag
            case 'individual'
                weight_vec = X(1:num_individuals); % just read input (how many null alleles in each individual)
            case 'summary' % We shouldn't even get here!
                weight_vec = []; % we don't have individual information
        end
    end
    for i=1:num_individuals          % Compute phenotype part. Loop on individuals. Heaviest loop
        log_like = log_like + ...
            log(internal_phenotype_fun(y(i), beta*weight_vec(i), ...
            sigma_e, trait_struct));
    end
end % switch if we have mixture or not



% Internal function: return the likelihood of phenotype given genotype
% where we know which alleles are functional. We use a Gaussian model (or
% liability-threshold for disease).
%
% Input:
% y - vector of phenotypes
% sum_x - total additive effect of all functional rare alleles for each person
% sigma_e - evnironmental noise
% trait_struct - type - quantitative or binary (disease)
%                prevalence - when trait is binary (disease)
%
% Output:
% ret - likelihood of phenotype for each individual
%
function ret = internal_phenotype_fun(y, sum_x, sigma_e, trait_struct)

if(~exist('trait_struct', 'var'))
    trait_struct = [];
    trait_struct.type = 'quantitative';
end

switch trait_struct.type
    case {'quantitative', 'continuous', 'QTL'}
        ret = normpdf( (y - sum_x) ./ sigma_e ); % For a gaussian quantitative trait
    case {'discrete', 'disease'}
        x_mu = norminv(1-trait_struct.prevalence);
        %        ret = 1-normcdf( (x_mu - sum_x) / sigma_e );
        ret = y + (1-2.*y) .* normcdf( (x_mu - sum_x) / sigma_e );  % use liability transformation for disease trait. take cumulative (tail) probability
end

% Inernal function for setting parameters 
function params = internal_set_parmas(params)


if(~isfield(params, 'include_phenotype') || isempty(params.include_phenotype))
    params.include_phenotype = 1; % flag saying if to include phenotypes when computing likelihood (default is one.)
end
if(params.include_phenotype == -1) % set if to copmute genotype part of likelihood
    params.include_genotype = 0;
else
    params.include_genotype = 1;
end

if(~isfield(params, 'null_w_vec') || isempty(params.null_w_vec))  % Perform exponential search on all values of w
    params.optimize_alpha = 1;
else  % here we know which alleles are null and which are neutral
    params.optimize_alpha = 0;
end
if(size(params.null_w_vec, 2) == 2) % get also true null state (for each missense allele)
    params.null_w_vec = params.null_w_vec(:,1);
end


if(~isfield(params, 'full_flag') || isempty(params.full_flag)) % default: use entire genotypes and phenotypes (not just summary statistics)
    params.full_flag = 1;
end

if(~isfield(params, 'print_flag') || isempty(params.print_flag)) % default: don't print
    params.print_flag = 0;
end
