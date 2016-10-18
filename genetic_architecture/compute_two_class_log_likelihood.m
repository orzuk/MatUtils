% Compute log-likelihood for data using the two-class model
% Formula is taken from eq. (XXX) from document on RVAS: 
%
% Input:
% s_null_vec - vector of possible selection coefficients for the null alleles. Should be NEGATIVE for delterious alleles (so fitness is 1+s)
% alpha_vec - vector of fraction of null alleles at birth
% beta_vec - vector of effect size of null alleles
% rare_cumulative_per_gene - total expected number of rare alleles in gene (related to theta)
% target_size_by_class_vec - NEW! target size for each allele type (used in poisson model). 
%                            This is like the above but counts expected number of different alleles, 
%                            not their frequencies !  %%% NEW! alleles class type vector. This states for each allele if it is neutral, harmfull, or in a mixture (missense)
% N - effective population size
% X - genotype data matrix - for some models we actually need it. For others, just the allele frequencies.
% y - phenotype data vector - (optional) may be used if likelihood includes phenotype 
% trait_type - disease or quantitative
% prevalence - frequency of disease in the population for disease trait
% null_w_vec - (optional) this is the assignment of which alleles are null and which not.
%                   This simplifies the likelihood computation a lot when it is known
%                   (no need to average over possible assignments of w). Convention:
%                      1: null
%                      0: neutral
%                     -1: missense (unknown, drawn from a mixture)
% include_phenotype - flag saying which part of the likelihood we compute:
%                      1: all (genotype+phenotype) (default)
%                      0: only genotype part of LL
%                     -1: only phenotype part of LL
% full_flag - input format:
%                      1: X is genotype matrix (default)
%                      0: X contains sufficient statistics (sums of rows and columns)
% num_individuals - input number of individuals (when not given in input vectors)
% D - demographic parameters (NEW! optional! (currently only constant pop. size supported)
%
% Output:
% log_like_mat - Matrix (3-d) of log-likelihood of data for each parameter choice (s, alpha and beta)
% P_poly - ??? Probability of polymorphic alleles? 
% 
function [log_like_mat, P_poly] = ...
    compute_two_class_log_likelihood(s_null_vec, alpha_vec, beta_vec, rare_cumulative_per_gene, target_size_by_class_vec, N, ...
    X, y, trait_type, prevalence, null_w_vec, include_phenotype, full_flag, num_individuals, D, print_flag)


%%%%%%%%%%%%%%%%%%%%%%%% Set Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
ttt=cputime; 
use_allele_freq_flag=2; % should be 2. (temp for debugging - allow computing only partial likelihoods) 
AssignGeneralConstants;
if(~exist('trait_type', 'var') || isempty(trait_type))
    trait_type = 'quantitative';
end
if(~exist('prevalence', 'var') || isempty(prevalence))
    prevalence = [];
end
if(~exist('include_phenotype', 'var') || isempty(include_phenotype))
    include_phenotype = 1; % flag saying if to include phenotypes when computing likelihood (default is one.)
end
if(include_phenotype == -1) % set if to copmute genotype part of likelihood
    include_genotype = 0;
else
    include_genotype = 1;
end
if(~include_phenotype) % here beta doesn't influence result
    beta_vec = 1;
end
if(~exist('target_size_by_class_vec', 'var') || isempty(target_size_by_class_vec))
    poisson_model_flag = 0; expand_format_flag = 'individual'; % include individual-level information (how many alleles)
else
    poisson_model_flag = 1; expand_format_flag = 'summary'; % use poisson model for each class
    target_size_neutral_alleles = target_size_by_class_vec(1); % convention: [neutral, null, missense]
    target_size_null_alleles = target_size_by_class_vec(2);
    target_size_missense_alleles = target_size_by_class_vec(3);
end
if(~exist('null_w_vec', 'var') || isempty(null_w_vec))                 % Perform exponential search on all values of w
    optimize_alpha = 1;
else  % here we know which alleles are null and which are neutral
    optimize_alpha = 0;
end
if(~exist('full_flag', 'var') || isempty(full_flag)) % default: use entire genotypes and phenotypes (not just summary statistics)
    full_flag = 1;
end
if(size(null_w_vec, 2) == 2) % get also true null state (for each missense allele) 
    true_null_w_vec = null_w_vec(:,2); 
    null_w_vec = null_w_vec(:,1); 
end
if(~exist('print_flag', 'var') || uisempty(print_flag)) % default: don't print 
    print_flag = 0; 
end
%%%%%%%%%%%%%%%%%%%%%%%% End Set Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


num_s = length(s_null_vec);
num_alpha = length(alpha_vec);
num_beta = length(beta_vec);
if(poisson_model_flag) % Compute counts
%    num_distinct_null_alleles_observed = sum(null_w_vec == 1);  % Polymorphic in POPULATION. these may be not polymorphic in sample
%    num_distinct_neutral_alleles_observed = sum(null_w_vec == 0); % Polymorphic in POPULATION. these may be not polymorphic in sample
%    num_distinct_missense_alleles_observed = sum(null_w_vec == -1); % mixture. % Polymorphic in POPULATION. these may be not polymorphic in sample
    num_polymorphic_null_alleles_observed = sum(X(null_w_vec == 1,:)>0); % polymorphic in SAMPLE
    num_polymorphic_neutral_alleles_observed = sum(X(null_w_vec == 0,:)>0); % polymorphic in SAMPLE
    num_polymorphic_missense_alleles_observed = sum(X(null_w_vec == -1,:)>0); % polymorphic in SAMPLE
%    if(exist('true_null_w_vec', 'var'))
%        num_polymorphic_missense_null_alleles_observed = sum(X((null_w_vec == -1)&(true_null_w_vec==1),:)>0); 
%        num_polymorphic_missense_neutral_alleles_observed = ...
%            num_polymorphic_missense_alleles_observed - num_polymorphic_missense_null_alleles_observed;
%    end
%    num_polymorphic_alleles_observed = num_polymorphic_null_alleles_observed + ...
%        num_polymorphic_neutral_alleles_observed+num_polymorphic_missense_alleles_observed; % total # of poymorphic alleles in sample (should replace L)
    if(~isfield(D, 'mu'))
        D.mu = 1.5 * 10^(-8); % TEMP! Set a default value for mu !
    end
    
    prob_null_allele_polymorphic_in_population = zeros(num_s, 1);
    prob_null_allele_polymorphic_in_sample = zeros(num_s, 1);
    prob_missense_allele_polymorphic_in_population = zeros(num_s, num_alpha);
    prob_missense_allele_polymorphic_in_sample = zeros(num_s, num_alpha);
end

% theta = rare_cumulative_per_gene; % total rate of polymorphic alleles. So far not used !!! 
sigma_e = 1; % environmental noise level
if(full_flag)
    [num_individuals, L] = size(X); % set number of individuals and number of SNPs
else
    if(~isempty(y)) % when phenotype is given
        num_individuals = length(y);
    else % here we have no way of knowing the # of individuals !!!
    end
    %    L = length(X) - num_individuals; % compute L later !!! set numbers from compact representation
end

log_like_mat = zeros(num_s, num_alpha, num_beta);
x_vec = (1:2*N-1) ./ (2*N); % vector of allele frequencies
if(full_flag)
    x_inds = cell(num_individuals,1);
    for i=1:num_individuals % record alternative alleles
        x_inds{i} = find(X(i,:));
    end
    % [unique_num_people unique_inds num_dup] = unique_with_inds(num_carriers_vec'); % why do we need this?
    num_carriers_vec = sum(X); % num. of individuals carrying each rare allele
    num_alleles_vec = sum(X,2); % num. of rare alleles in each individual
    num_individuals_vec = repmat(num_individuals, L, 1); % just repeat (assume for all alleles the same # of individuals is profiled)
else % here we only have summary statistics (sum of rows and columns)
    [num_alleles_vec, num_null_alleles_vec, num_carriers_vec, num_individuals_vec, L] = ...
        expand_two_class_summary_statistics(X, num_individuals, expand_format_flag);% get num. of individuals for each rare allele, num. of rare allele for each individual
    %     num_carriers_vec = X(end-L+1:end,:); % num. of individuals for each rare allele
    %     num_alleles_vec = X(1:num_individuals,:);  % num. of rare alleles in each individual
end
allele_freq_vec = num_carriers_vec ./ num_individuals; % get observed allele frequency of each derived allele

[unique_num_alleles_vec, I_num_alleles, J_num_alleles] = unique(num_alleles_vec);

neutral_allele_freq_hist = exp(allele_freq_spectrum(x_vec, 0, N, 0, 'log')); % allele freq. distribution for neutral alleles. NOT Normalized!
sum_neutral_allele_freq_hist = sum(neutral_allele_freq_hist);

prob_null_given_x = zeros(L,1); % conditional probability of allele being null when we know th
%%% T_0 = absorption_time_by_selection(0, 1, N, 1/(2*N), 1-1/(2*N), 0); % 'freq');
%%% T_0 = integral_hist(x_vec, neutral_allele_freq_hist);
%T_0 = absorption_time_by_selection(0, 1, N, 1/(2*N), 1-1/(2*N), 0); % use analytic approximation (not histogtam). Turns out to matter a lot!
T_0 = sum(neutral_allele_freq_hist) * (x_vec(2)-x_vec(1)); 

log_x_vec = log(x_vec); log_one_minus_x_vec = log(1-x_vec);
underflow_correction_p = min(num_individuals-1, max(1, num_carriers_vec)) ./ num_individuals;
underflow_correction_q = 1-underflow_correction_p;
log_like_correction = num_carriers_vec .* log(underflow_correction_p) + (num_individuals-num_carriers_vec) .* log(underflow_correction_q); % update log-likelihood: correction for multinomial/binomial sampling?
tmp_likelihood_one_allele = repmat(BIG_NUM, max(num_individuals_vec), 3);

if(optimize_alpha) % Here what ?? 
    cond_y_tab = zeros(num_individuals, max(num_alleles_vec)+1); % table. entry (i,j) is Prob. (y(i) | beta*(j-1) nulls)
    if(include_phenotype) % prepare binomial tables
        binom_vec = zeros(length(unique_num_alleles_vec), max(unique_num_alleles_vec)+1); % save binomial coefficients in table
    end
end
if(poisson_model_flag)
    lambda_s = zeros(num_s, 1);
    lambda_missense = zeros(num_s, num_alpha);
    tmp_z0_vec = exp( num_individuals_vec(1) .* (log_one_minus_x_vec)); %  - log(underflow_correction_q(1))  ); % what's this? multinomial/binomial coefficient?
    prob_neutral_allele_polymorphic_in_population = 4 * N * D.mu * T_0; % prob. allele polymorphic in population
    prob_neutral_allele_polymorphic_in_sample = prob_neutral_allele_polymorphic_in_population .* ...
        (1 - sum( neutral_allele_freq_hist .* (1-x_vec) .^ num_individuals_vec(1) ) / sum(neutral_allele_freq_hist));
    % (1 - integral_hist(x_vec,  neutral_allele_freq_hist .* tmp_z0_vec ) ./ integral_hist(x_vec, neutral_allele_freq_hist) ); % This gives wrong results!!
    
    [unique_num_carriers, I, J] = unique(num_carriers_vec);
    tmp_z_vec = sparse(length(unique_num_carriers), length(x_vec)); % creat sparse matrix
    pos_tmp_z_inds = cell(length(unique_num_carriers), 1);
    %%    t2 = cputime; % tmp_ind_vec = []; tmp_val_vec = [];
    for j=1:length(unique_num_carriers) % loop on unique 
%        run_j = j 
        tmp_vec = ...
            exp( log_binom(num_individuals_vec(I(j)), num_carriers_vec(I(j))) + ...
            num_carriers_vec(I(j)) .* log_x_vec + ...
            (num_individuals_vec(I(j))-num_carriers_vec(I(j))) .* log_one_minus_x_vec ); % assumes num individuals is the same
        pos_tmp_z_inds{j} = find(tmp_vec > 10^(-8));
        tmp_z_vec(j,pos_tmp_z_inds{j}) = tmp_vec(pos_tmp_z_inds{j});
    end
    %%    cputime - t2
    
end


T_s = zeros(num_s,1); % time an allele spends at polymorphic state 
for i_s = 1:num_s % loop on parameters

    % Get allele frequency based on spectrum 
    null_allele_freq_hist = exp(allele_freq_spectrum(x_vec, s_null_vec(i_s), N, 0, 'log')); % allele freq. distribution for null alleles.  NOT Normalized!
    sum_null_allele_freq_hist = sum(null_allele_freq_hist);
    %    T_s = integral_hist(x_vec, null_allele_freq_hist); %%% T_s = absorption_time_by_selection(-s_null_vec(i_s), 1, N, 1/(2*N), 1-1/(2*N), 0); % 'freq'); % sure we need to use 'freq' here ???
    %%    T_s(i_s) = absorption_time_by_selection(-s_null_vec(i_s), 1, N, 1/(2*N), 1-1/(2*N), 0); % use analytic approximation (not histogtam). Turns out to matter a lot!
    T_s(i_s) = sum(null_allele_freq_hist) * (x_vec(2)-x_vec(1)); 
    
    if(poisson_model_flag) % compute poisson part of likelihood
        % we must assume here we know the number of individuals profiled at each region
        %        tmp_z0_vec = exp( num_individuals_vec(1) .* (log_one_minus_x_vec)); %  - log(underflow_correction_q(1))  ); % what's this? multinomial/binomial coefficient?
        prob_null_allele_polymorphic_in_population(i_s) = 4 * N * D.mu * T_s(i_s); % prob. allele polymorphic in population
        prob_null_allele_polymorphic_in_sample(i_s) = prob_null_allele_polymorphic_in_population(i_s) .* ...
            (1 - sum( null_allele_freq_hist .* (1-x_vec) .^ num_individuals_vec(1) ) / sum(null_allele_freq_hist));
        %            (1 - integral_hist(x_vec,  null_allele_freq_hist .* tmp_z0_vec ) ./ integral_hist(x_vec, null_allele_freq_hist) );
        prob_missense_allele_polymorphic_in_population(i_s,:) = alpha_vec .* prob_null_allele_polymorphic_in_population(i_s) + ...
            (1-alpha_vec) .* prob_neutral_allele_polymorphic_in_population;
        prob_missense_allele_polymorphic_in_sample(i_s,:) = alpha_vec .* prob_null_allele_polymorphic_in_sample(i_s) + ...
            (1-alpha_vec) .* prob_neutral_allele_polymorphic_in_sample;
        
        % set lambdas for poisson model
        lambda_s(i_s) = 4 * N * D.mu * target_size_null_alleles * T_s(i_s);
        %        lambda_0 = 4 * N * D.mu * target_size_neutral_alleles * T_0;
        lambda_missense(i_s,:) = 4 * N * D.mu * target_size_missense_alleles .* ( alpha_vec .* T_s(i_s) + (1-alpha_vec) .* T_0 );
        
    end
% % % % %     p_null_vec_in_population = alpha_vec .* T_s(i_s) ./ (alpha_vec .* T_s(i_s) + (1-alpha_vec) .* T_0); % compute probability that a given observed polymorphic locus is null in the population!!!!
% % % % %     p_null_vec_in_sample = alpha_vec .* prob_null_allele_polymorphic_in_sample(i_s) ./ ...
% % % % %         (alpha_vec .* prob_null_allele_polymorphic_in_sample(i_s) + ...
% % % % %         (1-alpha_vec) .* prob_neutral_allele_polymorphic_in_sample);    
% % % % %     p_null_empiric = num_polymorphic_missense_null_alleles_observed / num_polymorphic_missense_alleles_observed;
% % % % %     figure; plot(alpha_vec, p_null_vec_in_population, '.'); hold on;
% % % % %     plot(alpha_vec, p_null_vec_in_sample, 'g.');
% % % % %     line([0 1], [p_null_empiric p_null_empiric], 'color', 'r');
% % % % %     xlabel('\alpha'); ylabel('Pr(null)'); legend('population', 'sample', 'empiric');
    for i_alpha = 1:num_alpha
        ttt = cputime;
        %        p_null = alpha_vec(i_alpha) * T_s(i_s) / (alpha_vec(i_alpha) * T_s(i_s) + (1-alpha_vec(i_alpha)) * T_0); % compute probability that a given observed polymorphic locus is null
        missense_allele_freq_hist = alpha_vec(i_alpha) .* null_allele_freq_hist + (1-alpha_vec(i_alpha)) .* neutral_allele_freq_hist; % mixture of neutral and null freq. distributions.  NOT Normalized! the weight of neutral and null may be very different!
        sum_missense_allele_freq_hist = sum(missense_allele_freq_hist);
        if(include_phenotype && optimize_alpha) % prepare binomial tables
            for i=1:length(unique_num_alleles_vec)
                binom_vec(i,1:unique_num_alleles_vec(i)+1) = ...
                    binopdf(0:unique_num_alleles_vec(i), unique_num_alleles_vec(i), alpha_vec(i_alpha));
            end
        end % if include phenotype
        
                        
        if(poisson_model_flag) % compute poisson part of likelihood  (observation of polymorphic alleles)           
            if( ismember(use_allele_freq_flag, [0, 2]))
                log_like_mat(i_s,i_alpha,:) = ...
                    (target_size_null_alleles - num_polymorphic_null_alleles_observed) .* ...
                    log(1-prob_null_allele_polymorphic_in_sample(i_s)) + ...
                    (target_size_neutral_alleles - num_polymorphic_neutral_alleles_observed) .* ... % Compute log-likelihood for alleles not present
                    log(1-prob_neutral_allele_polymorphic_in_sample) + ...
                    (target_size_missense_alleles - num_polymorphic_missense_alleles_observed) .* ...
                    log(1-prob_missense_allele_polymorphic_in_sample(i_s, i_alpha));
                
                log_like_mat(i_s,i_alpha,:) = log_like_mat(i_s,i_alpha,:) + ... % add log-likelihood for observed alleles 
                    num_polymorphic_null_alleles_observed .* ...
                    log(prob_null_allele_polymorphic_in_sample(i_s)) + ...
                    num_polymorphic_neutral_alleles_observed .* ...
                    log(prob_neutral_allele_polymorphic_in_sample) + ...
                    num_polymorphic_missense_alleles_observed .* ...
                    log(prob_missense_allele_polymorphic_in_sample(i_s, i_alpha));            % all probs here should be in POPULATION!!! (if we include all parts of likelihood!)
                
                
                if(use_allele_freq_flag == 0)
                    continue;
                end
            else
                log_like_mat(i_s,i_alpha,:) = 0;
            end            
            
            if( ismember(use_allele_freq_flag, [1, 2])) % add log of binomial coefficient (?)
                log_like_mat(i_s,i_alpha,:) = log_like_mat(i_s,i_alpha,:) + ...
                    num_polymorphic_null_alleles_observed .* ...
                    log(prob_null_allele_polymorphic_in_population(i_s)) + ...
                    num_polymorphic_neutral_alleles_observed .* ...
                    log(prob_neutral_allele_polymorphic_in_population) + ...
                    num_polymorphic_missense_alleles_observed .* ...
                    log(prob_missense_allele_polymorphic_in_population(i_s, i_alpha));            % all probs here should be in POPULATION!!! (if we include all parts of likelihood!)
                
                log_like_mat(i_s,i_alpha,:) = log_like_mat(i_s,i_alpha,:) - ...
                    num_polymorphic_null_alleles_observed .* ...
                    log(prob_null_allele_polymorphic_in_sample(i_s)) - ...
                    num_polymorphic_neutral_alleles_observed .* ...
                    log(prob_neutral_allele_polymorphic_in_sample) - ...
                    num_polymorphic_missense_alleles_observed .* ...
                    log(prob_missense_allele_polymorphic_in_sample(i_s, i_alpha));            % all probs here should be in POPULATION!!! (if we include all parts of likelihood!)
%                                 log_like_mat(i_s,i_alpha,:) = 0;
            end % if use allele freq flag
            
        end % if poisson model flag 
        for i_beta = 1:num_beta % loop on effect size
            tmp_likelihood_one_allele(:) = BIG_NUM; % repmat(BIG_NUM, max(num_individuals_vec), 3);
            %            time_before_looping_on_alleles = cputime-ttt
            
            for j=1:L    % Compute genotype part. Loop on loci. Compute for each # of carriers only once 
                if( (num_carriers_vec(j) == 0) || (null_w_vec(j) == 0) )% loop only on poylmorphic non-synonmous alleles
                    continue;
                end
                    
                num_carriers_ind = find(num_carriers_vec(j) == unique_num_carriers);
                if(include_genotype || optimize_alpha)
                    if(tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2 ) == BIG_NUM) % this means we've already computed for this frequency
                        
                        %                    if(~poisson_model_flag) % old version (not clear what is computed here)
                        % % % %                     tmp_z_vec = exp( num_carriers_vec(j) .* (log_x_vec  - log(underflow_correction_p(j)) ) + ...
                        % % % %                         (num_individuals_vec(j)-num_carriers_vec(j)) .* (log_one_minus_x_vec - log(underflow_correction_q(j)) ) ); % what's this? multinomial/binomial coefficient?
                        % % % %                     if(any(isnan(tmp_z_vec) | isinf(tmp_z_vec)))
                        % % % %                         problem_with_nans = 9999999999
                        % % % %                         %                        return;
                        % % % %                     end
                        
                        if(include_genotype)
                            switch null_w_vec(j)
                                case 1 % null (depends on s, not on alpha)
                                    %                             fff  = find(cur_allele_freq_hist);
                                    tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2) = ...
                                        log( max(10^(-100), sum(null_allele_freq_hist(pos_tmp_z_inds{num_carriers_ind}) .* ...
                                        tmp_z_vec(num_carriers_ind,pos_tmp_z_inds{num_carriers_ind})) / sum_null_allele_freq_hist) );
                                    %                                    cur_allele_freq_hist = null_allele_freq_hist; sum_cur_allele_freq_hist = sum_null_allele_freq_hist;
                                case 0 % neutral (doesn't depend on either s or alpha)
                                    tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2) = ...
                                        log( max(10^(-100), sum(neutral_allele_freq_hist(pos_tmp_z_inds{num_carriers_ind}) .* ...
                                        tmp_z_vec(num_carriers_ind,pos_tmp_z_inds{num_carriers_ind})) / sum_neutral_allele_freq_hist) );
                                    %                                    cur_allele_freq_hist = neutral_allele_freq_hist; sum_cur_allele_freq_hist = sum_neutral_allele_freq_hist;
                                case -1 % missense (depends on both s and alpha)
                                    tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2) = ...
                                        log( max(10^(-100), sum(missense_allele_freq_hist(pos_tmp_z_inds{num_carriers_ind}) .* ...
                                        tmp_z_vec(num_carriers_ind,pos_tmp_z_inds{num_carriers_ind})) / sum_missense_allele_freq_hist) );
                                    %                                    cur_allele_freq_hist = missense_allele_freq_hist; sum_cur_allele_freq_hist = sum_missense_allele_freq_hist;
                            end % new! allow different alleles
                            
                            %                             tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2) = ...
                            %                                 log( max(10^(-100), sum(cur_allele_freq_hist .* ...
                            %                                 (1-x_vec).^ (num_individuals_vec(j)-num_carriers_vec(j)) .* ...
                            %                                 x_vec.^  num_carriers_vec(j)) / sum(cur_allele_freq_hist)) );  % very slow !!!
                            %                             exp( num_carriers_vec(j) .* log_x_vec + ...
                            %                                 (num_individuals_vec(j)-num_carriers_vec(j)) .* log_one_minus_x_vec )) / sum_cur_allele_freq_hist) );
                            
                            
                            %%%%%%%%                            log( max(10^(-100), sum(cur_allele_freq_hist .* tmp_z_vec) / (2*N)) ); % sum(cur_allele_freq_hist)) );
                            %%%%%%%%                        integral_hist(x_vec,  cur_allele_freq_hist .* tmp_z_vec )) ); % + ...
                            
                            
                            %                    num_carriers_vec(j) * log(underflow_correction_p(j)) + (n-num_carriers_vec(j)) * log(underflow_correction_q(j)); % update log-likelihood
                            %                 log_like_vec{i_s}(j) = log( integral_hist(x_vec,  (alpha_vec(i_alpha) .* null_allele_freq_hist + (1-alpha_vec(i_alpha)) .* neutral_allele_freq_hist) .* ...
                            %                     tmp_z_vec ) ) + ...
                            %                     num_carriers_vec(j) * log(underflow_correction_p(j)) + (n-num_carriers_vec(j)) * log(underflow_correction_q(j)); % update log-likelihood
                        end
                    end % already computed likelihood
                    log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                        tmp_likelihood_one_allele(num_carriers_vec(j), null_w_vec(j)+2);
                    
                    
                end % include genotypes or optimize alpha
                if(mod(j, 1000) == 0)
                    run_locus = j
                end
                if(include_phenotype && optimize_alpha)  % compute conditional probability of allele being null given frequency x.
                    prob_null_given_x(j) = p_null * integral_hist(x_vec, null_allele_freq_hist .* tmp_z_vec);
                    prob_null_given_x(j) = prob_null_given_x(j) / (prob_null_given_x(j) + ...
                        (1-p_null) * integral_hist(x_vec, neutral_allele_freq_hist .* tmp_z_vec));
                end
            end % loop on loci
            %            ttt_loop_on_loci = cputime - ttt
            if(include_genotype) % here add genotype part of log-likelihood 
                log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + 0 * sum(log_like_correction); % add correction. Just a constant (who cares)
                
                if(poisson_model_flag) % here include all alleles classes
                    log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) - 0; %%%%... % add normalization constant
                    %%%                        num_distinct_neutral_alleles_observed * log(T_0) - num_distinct_null_alleles_observed * log(T_s(i_s)) - ...
                    %%%                        num_distinct_missense_alleles_observed * log( alpha_vec(i_alpha) * T_s(i_s) + (1-alpha_vec(i_alpha)) * T_0 );
                else % here assume all alleles are missnese
                    log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) - ... % add normalization constant
                        L * log( alpha_vec(i_alpha) * T_s(i_s) + (1-alpha_vec(i_alpha)) * T_0 );
                end
                
            end % if include genotypes 
            
            % Skip this for now (we just use genotypes)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% End if include phenotypes %%%%%%%%%%%%
            if(include_phenotype) % compute phenotypes contribution to likelihood
                switch trait_type
                    case {'binary', 'disease'}
                        mean_x = mean(X(:));
                        mean_f = mean_x * alpha_vec(i_alpha);
                        var_exp = beta_vec(i_beta).^2 * mean_f * (1-mean_f);
                        sigma_e = 1 - var_exp;
                end
                if(optimize_alpha)       % Perform exponential search on all values of w (which alleles are null and which not). Is this the right formula?
                    for i=1:n
                        %                        for j=1:num_alleles_vec(i)+1
                        cond_y_tab(i,1:num_alleles_vec(i)+1) = internal_phenotype_fun( ...
                            y(i), beta_vec(i_beta).*(0:num_alleles_vec(i)), sigma_e, trait_type, prevalence);
                        %                        end
                    end
                    if(full_flag) % exponential search
                        W_mat = my_dec2base( 0:2^max(num_alleles_vec)-1, 2, max(num_alleles_vec)); % The heavy part: exponentially many genotypes
                        W_sum_cell = cell(max(num_alleles_vec)+1,1); % what is this?
                        for i=0:max(num_alleles_vec)
                            W_sum_cell{i+1} = sum(W_mat(1:2^i,1:i),2);
                        end
                        
                        for i=1:n          % Compute phenotype part. Loop on individuals. Heaviest loop
                            %                 i_is = i
                            %                 num_alleles_is =num_alleles_vec(i)
                            k = num_alleles_vec(i);
                            prob_null_mat = repmat(prob_null_given_x(x_inds{i}), 1, 2^k)';
                            prob_null_mat = prod( prob_null_mat .^ W_mat(1:2^k,1:k) .* (1-prob_null_mat) .^ (1-W_mat(1:2^k,1:k)), 2);
                            log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                                log(sum(prob_null_mat .* cond_y_tab(i, W_sum_cell{k+1}+1)')); % why plus one???
                        end % loop on individuals
                    else % use binomial distribution (linear search instead of exponential) % here X is summary statistics (row and column sum)
                        for i=1:n % loop on individuals (heavy part)
                            log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                                log( sum( binopdf(0:num_alleles_vec(i), num_alleles_vec(i), alpha_vec(i_alpha)) .* ...
                                cond_y_tab(i,1:num_alleles_vec(i)+1)  ) ); % sum over binomial probabilities
                            log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                                log( sum( binom_vec(J_num_alleles(i), 1:num_alleles_vec(i)+1) .* ...
                                cond_y_tab(i,1:num_alleles_vec(i)+1)  ) ); % sum over binomial probabilities
                            
                        end
                        
                    end
                else % here assume that we know which alleles are null and which are not
                    if(full_flag) % here X is a genotype matrix
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
                        log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                            log(internal_phenotype_fun(y(i), beta_vec(i_beta)*weight_vec(i), ...
                            sigma_e, trait_type, prevalence));
                    end
                end % switch if we have mixture or not
            end % if to include phenotypes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% End if include phenotypes %%%%%%%%%%%%
            if(print_flag)
                if( (mod(i_alpha, 50) == 0) || mod(i_s, 50) == 0)
                    run_index_s_alpha_beta = [i_s i_alpha i_beta]
                    time_one_likelihood = cputime-ttt
                end
            end
        end % loop on effect size beta
        
            
    end % loop on mixture coefficient alpha
end % loop on selection coefficients

P_poly = var2struct(prob_neutral_allele_polymorphic_in_population, prob_neutral_allele_polymorphic_in_sample, ...
    prob_null_allele_polymorphic_in_population, prob_null_allele_polymorphic_in_sample, ...
    prob_missense_allele_polymorphic_in_population, prob_missense_allele_polymorphic_in_sample);



% Internal function: return the likelihood of phenotype given genotype
% where we know which alleles are functional. We use a Gaussian model (or
% liability-threshold for disease).
%
% Input:
% y - vector of phenotypes
% sum_x - total additive effect of all functional rare alleles for each person
% sigma_e - evnironmental noise
% trait_type - quantitative or binary (disease)
% prevalence - when trait is binary (disease)
%
% Output:
% ret - likelihood of phenotype for each individual
%
function ret = internal_phenotype_fun(y, sum_x, sigma_e, trait_type, prevalence)

if(~exist('trait_type', 'var'))
    trait_type = 'quantitative';
end

switch trait_type
    case {'quantitative', 'continuous', 'QTL'}
        ret = normpdf( (y - sum_x) ./ sigma_e ); % For a gaussian quantitative trait
    case {'discrete', 'disease'}
        x_mu = norminv(1-prevalence);
        %        ret = 1-normcdf( (x_mu - sum_x) / sigma_e );
        ret = y + (1-2.*y) .* normcdf( (x_mu - sum_x) / sigma_e );  % use liability transformation for disease trait. take cumulative (tail) probability
end


