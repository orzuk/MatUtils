% Compute log-likelihood for data using the two class model
%
% Input:
% s_null_vec - vector of possible selection coefficients for the null alleles
% alpha_vec - vector of fraction of null alleles at birth
% beta_vec - vector of effect size of null alleles
% rare_cumulative_per_gene - total expected number of rare alleles in gene (related to theta)
% N - effective population size
% X - genotype data matrix
% y - phenotype data vector
% trait_type - disease or quantitative
% prevalence - frequencu of disease for disease trait
% null_w_vec - (optional) this is the assignment of which alleles are null
%            and which not. This simplifies the likelihood computation alot when it is known
%
% Output:
% log_like_mat - Matrix (3-d) of log-likelihood of data for each parameter choice
%
function log_like_mat = ...
    compute_two_class_likelihood(s_null_vec, alpha_vec, beta_vec, rare_cumulative_per_gene, N, ...
    X, y, trait_type, prevalence, null_w_vec)


if(~exist('trait_type', 'var') || isempty(trait_type))
    trait_type = 'quantitative';
end
if(~exist('prevalence', 'var') || isempty(prevalence))
    prevalence = [];
end

include_phenotype = 1; % flag saying if to include phenotypes when computing likelihood (default is one.)

theta = rare_cumulative_per_gene;
sigma_e = 1; % environmental noise
[n L] = size(X); % set number of individuals and number of SNPs

num_s = length(s_null_vec);
num_alpha = length(alpha_vec);
num_beta = length(beta_vec);
log_like_mat = zeros(num_s, num_alpha, num_beta);
x_vec = (1:2*N-1) ./ (2*N);

num_people_vec = sum(X); % num. of individuals for each rare allele
num_alleles_vec = sum(X,2); % num. of rare alleles in each individual

x_neutral_hist = exp(allele_freq_spectrum(x_vec, 0, N, 0, 'log'));
prob_null_given_x = zeros(L,1); % conditional probability of allele being null
t_0 = absorbtion_time_by_selection(0, 1, N, 1/(2*N), 1-1/(2*N), 0); % 'freq');
p_neutral_vec = exp(allele_freq_spectrum(x_vec, 0, N, 0, 'log')); % allele freq. distribution for neutral alleles

x_inds = cell(n,1);
for i=1:n
    x_inds{i} = find(X(i,:));
end
[unique_num_people unique_inds num_dup] = unique_with_inds(num_people_vec');

log_x_vec = log(x_vec); log_one_minus_x_vec = log(1-x_vec);
underflow_correction_p = min(n-1, max(1, num_people_vec)) ./ n;
underflow_correction_q = 1-underflow_correction_p;
log_like_correction = num_people_vec .* log(underflow_correction_p) + (n-num_people_vec) .* log(underflow_correction_q); % update log-likelihood
for i_s = 1:num_s
    p_null_vec = exp(allele_freq_spectrum(x_vec, s_null_vec(i_s), N, 0, 'log')); % allele freq. distribution for null alleles
    t_s = absorbtion_time_by_selection(-s_null_vec(i_s), 1, N, 1/(2*N), 1-1/(2*N), 0); % 'freq'); % sure we need to use 'freq' here ???
    x_null_hist = exp(allele_freq_spectrum(x_vec, s_null_vec(i_s), N, 0, 'log'));
    
    for i_alpha = 1:num_alpha
        p_null = alpha_vec(i_alpha) * t_s / (alpha_vec(i_alpha) * t_s + (1-alpha_vec(i_alpha)) * t_0); % compute probability that a given locus is null
        tmp_z_mix_vec = alpha_vec(i_alpha) .* x_null_hist + (1-alpha_vec(i_alpha)) .* x_neutral_hist;
        for i_beta = 1:num_beta
            run_index = [i_s i_alpha i_beta]
            %           log_like_mat(i_s,i_alpha,i_beta)=0;
            for j=1:L          % Compute genotype part. Loop on loci
                tmp_z_vec = exp( num_people_vec(j) .* (log_x_vec  - log(underflow_correction_p(j)) ) + ...
                    (n-num_people_vec(j)) .* (log_one_minus_x_vec - log(underflow_correction_q(j)) ) );
                log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                    log( integral_hist(x_vec,  tmp_z_mix_vec .* tmp_z_vec ) ); % + ...
                %                    num_people_vec(j) * log(underflow_correction_p(j)) + (n-num_people_vec(j)) * log(underflow_correction_q(j)); % update log-likelihood
                %                 log_like_vec{i_s}(j) = log( integral_hist(x_vec,  (alpha_vec(i_alpha) .* x_null_hist + (1-alpha_vec(i_alpha)) .* x_neutral_hist) .* ...
                %                     tmp_z_vec ) ) + ...
                %                     num_people_vec(j) * log(underflow_correction_p(j)) + (n-num_people_vec(j)) * log(underflow_correction_q(j)); % update log-likelihood
                if(mod(j, 100) == 0)
                    j_is = j
                end
                if(include_phenotype)
                    prob_null_given_x(j) = p_null * integral_hist(x_vec, p_null_vec .* tmp_z_vec);
                    prob_null_given_x(j) = prob_null_given_x(j) / (prob_null_given_x(j) + ...
                        (1-p_null) * integral_hist(x_vec, p_neutral_vec .* tmp_z_vec));
                end
                %                0.5; % compute conditional probability of allele being null given x. Should change this !!!!
            end % loop on loci
            log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + sum(log_like_correction); % add correction
            log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) - ... % add normalization constant
                L * log( alpha_vec(i_alpha) * absorbtion_time_by_selection(-s_null_vec(i_s), theta, N, x_vec(1), x_vec(end), 0) + ...
                (1-alpha_vec(i_alpha)) * absorbtion_time_by_selection(0, theta, N, x_vec(1), x_vec(end), 0) );
            
            
            if(include_phenotype) % here the phenotypes are also considered for disease
                switch trait_type
                    case {'binary', 'disease'}
                        mean_x = mean(X(:));
                        mean_f = mean_x * alpha_vec(i_alpha);
                        var_exp = beta_vec(i_beta).^2 * mean_f * (1-mean_f);
                        sigma_e = 1 - var_exp;
                end
                
                cond_y_tab = zeros(n, max(num_alleles_vec)+1); % table. entry (i,j) is Prob. (y | beta*(j-1) nulls)
                
                % Perform exponential search on all values of w
                if(~exist('null_w_vec', 'var') || isempty(null_w_vec))
                    for i=1:n
                        for j=1:num_alleles_vec(i)+1
                            cond_y_tab(i,j) = internal_phenotype_fun(y(i), beta_vec(i_beta)*(j-1), ...
                                sigma_e, trait_type, prevalence);
                        end
                    end
                    W_mat = my_dec2base( 0:2^max(num_alleles_vec)-1, 2, max(num_alleles_vec)); % The heavy part:
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
                else % here assume that we know which alleles are null and which are not
                    weight_vec = X * vec2column(null_w_vec); % vector saying how many null alleles are in each individual
                    for i=1:n          % Compute phenotype part. Loop on individuals. Heaviest loop
                        log_like_mat(i_s,i_alpha,i_beta) = log_like_mat(i_s,i_alpha,i_beta) + ...
                            log(internal_phenotype_fun(y(i), beta_vec(i_beta)*weight_vec(i), ...
                            sigma_e, trait_type, prevalence));
                    end
                end
                
            end % temp include phenotypes
        end
    end
end % loop on seleection coefficients


% figure; plot(log_like_vec{1}, log_like_vec{2}, '.'); hold on;
% plot(min(log_like_vec{1}):max(log_like_vec{1}), min(log_like_vec{1}):max(log_like_vec{1}), 'r');
% xlabel(['s=' num2str(s_null_vec(1))]); ylabel(['s=' num2str(s_null_vec(2))]);
% XXX = 324


% Internal function: return the likelihood of phenotype given genotype
%
% Input:
% y - vector of phenotypes
% sum_x - total additive effect of all functional rare alleles for each person
% sigma_e - evniromnetal noise
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
    case {'quantitative', 'continuous'}
        ret = normpdf( (y - sum_x) ./ sigma_e ); % For a gaussian quantitative trait
    case {'discrete', 'disease'}
        x_mu = norminv(1-prevalence);
        ret = 1-normpdf( (x_mu - sum_x) / sigma_e );
end


