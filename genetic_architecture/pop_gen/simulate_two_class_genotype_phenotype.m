% Simulate rare alleles data according to selection coefficient
%
% Input:
% s_null - selection coefficient of null alleles (NEGATIVE for deleterious alleles)
% alpha - fraction of null alleles
% beta - effect size of null alleles on trait
% rare_cumulative_per_gene - total fraction of polymorphic alleles
% N - population size
% L - target size in nucleotides
% num_individuals - number of individuals
% iters - number of iterations to simulate
% trait_type - binary (disease) or continuous
% prevalence - for disease traits
% full_flag - saying if to output all individuals' genotypes (default) or just summary statistics
% poisson_flag - NEW! allow simulation of multiple classes !
%
% Output:
% X - Genotype data matrix
% y - Phenotype data vector
% is_null_mat - vector saying for each allele if it was null (harmful)
% f_mat - vector of derived allele frequencies
% allele_type_vec - NEW (optional) vector saying fro each allele if it is synonymous, missense, or stop
%
function [X, y, is_null_mat, f_mat, allele_type_vec, P_poly] = simulate_two_class_genotype_phenotype(s_null, alpha, beta, ...
    rare_cumulative_per_gene, target_size_by_class_vec, ...
    N, L, num_individuals, iters, trait_type, prevalence, full_flag, poisson_flag)

Assign24MammalsGlobalConstants; AssignRVASConstants;
mu= mu_per_site;

% sigma_e = 1; % noise due to enviromnent and other hidden genetic effects
if(~exist('full_flag', 'var') || isempty(full_flag)) % simulate entire genotype
    full_flag = 1;
end
if(~exist('poisson_flag', 'var') || isempty(poisson_flag))
    poisson_flag = 0;
end
t_s = absorption_time_by_selection(s_null, 1, N, 1/(2*N), 1-1/(2*N), 0);
t_0 = absorption_time_by_selection(0, 1, N, 1/(2*N), 1-1/(2*N), 0);
p_null = alpha * t_s / (alpha * t_s + (1-alpha) * t_0); % compute probability that each locus is null

if(poisson_flag) % randomize number of alleles
    lambda_by_class_vec = 4.*N.*mu .* [t_0 t_s (alpha*t_s + (1-alpha)*t_0)] .*  target_size_by_class_vec;  % format: [synonymous, stop, missense]
    num_polymorphic_alleles_by_class_vec = zeros(length(lambda_by_class_vec), iters);
    for i=1:iters  % randomize number of total polymorphic alleles
        num_polymorphic_alleles_by_class_vec(:,i) = vec2column(poissrnd(lambda_by_class_vec));
    end
    max_num_polymorphic_alleles_by_class_vec = max(num_polymorphic_alleles_by_class_vec,[],2); % get maximum number of alleles
    L = sum(max_num_polymorphic_alleles_by_class_vec); % Get total alleles
    expand_format_flag = 'summary';
else
    expand_format_flag = 'individual';
end
if(full_flag) % prepare genotype matrix and phenotype vector
    X = zeros(num_individuals, L, iters);
else
    switch expand_format_flag
        case 'summary'
            X = zeros(2*L, iters); % just keep sum to save space. Format: [sum-on-loci;  num-individuals-profiled]
        case 'individual'
            X = zeros(2*num_individuals+L, iters); % just keep sum of rows and columns to save space. Format: [sum-on-loci;  sum-null-on-loci; sum-on-individuals]
    end
end
y = zeros(num_individuals, iters);


if(poisson_flag) % set which alleles are null and which not
    is_null_mat = rand(L, iters) < p_null; % set for each polymorphic allele whether it is null or not
    allele_type_vec = zeros(L,iters);
    for i=1:iters
        is_null_mat(1:num_polymorphic_alleles_by_class_vec(1,i),i) = 0; % all neutral alleles
        is_null_mat((sum(num_polymorphic_alleles_by_class_vec(1,i))+1):sum(num_polymorphic_alleles_by_class_vec(1:2,i)),i) = 1; % all null alleles
        allele_type_vec(1:num_polymorphic_alleles_by_class_vec(1,i),i) = SYNONYMOUS;
        allele_type_vec((num_polymorphic_alleles_by_class_vec(1,i)+1):sum(num_polymorphic_alleles_by_class_vec(1:2,i)),i) = STOP;
        allele_type_vec((sum(num_polymorphic_alleles_by_class_vec(1:2,i))+1):sum(num_polymorphic_alleles_by_class_vec(:,i)),i) = MISSENSE;
    end % loop on iters
else
    is_null_mat = rand(L, iters) < p_null; % set for each allele whether it is null or not
    allele_type_vec = repmat(MISSENSE, L, 1); % all are missnese
end % if poisson
f_mat = zeros(L, iters); % Get allele frequencies

x_vec = (1:2*N-1)./(2*N);
p_null_vec = exp(allele_freq_spectrum(x_vec, s_null, N, 0, 'log'));
p_neutral_vec = exp(allele_freq_spectrum(x_vec, 0, N, 0, 'log'));


var_explained = zeros(iters,1);
for i=1:iters % simulate many times allele frequencies
    f_mat(:,i) = distrnd(x_vec, p_null_vec, L, 1) .* is_null_mat(:,i) + ...
        distrnd(x_vec, p_neutral_vec, L, 1) .* (1-is_null_mat(:,i));  % sample allele frequencies (need to add some computation here!!!)
    
    var_explained(i) = beta^2 * sum(f_mat(:,i) .* (1-f_mat(:,i)));  % Determine variance explained
    sigma_e = 1 - sqrt(var_explained(i));
    
    temp_x = rand(num_individuals, L) < repmat(f_mat(:,i), 1, num_individuals)'; % sample genotypes according to frequencies
    y(:,i) = sum( beta * temp_x .* repmat(is_null_mat(:,i), 1, num_individuals)', 2 ) + randn(num_individuals, 1) * sigma_e; % simulate phenotype conditional on genotype
    if(full_flag)
        X(:,:,i) = temp_x;
    else
        switch expand_format_flag
            case 'individual'
                X(1:num_individuals,i) = sum(temp_x,2); % num. of alleles in each individual
                X(num_individuals+1:2*num_individuals,i) = sum(temp_x .* repmat(is_null_mat(:,i)', num_individuals, 1),2); % num. of null alleles in each individual
                X(2*num_individuals+1:2*num_individuals+L,i) = sum(temp_x); % num. of individual carriers for each allele
            case 'summary'
                X(1:L,i) = sum(temp_x); % number of carriers
                X((L+1):2*L,i) = num_individuals; % currently assume number of individuals equal for all alleles
        end
    end % if full flag
end % loop on iters

if(~exist('trait_type', 'var'))
    trait_type = 'quantitative';
end
switch trait_type % transform to binary
    case {'disease', 'binary'}
        x_mu = norminv(1-prevalence);
        y = NormalizeData(y) > x_mu; % transform to binary r.v.s. (using threshold)
        %        y = ((y-mean(y))./std(y)) > x_mu; % transform to binary r.v.s.
end

if(poisson_flag)
    total_num_alleles = size(X, 1) / 2;    
    empirical_prob_neutral_polymorphic_in_population = sum(allele_type_vec(:)==SYNONYMOUS) / (target_size_by_class_vec(1) * iters);
    empirical_prob_null_polymorphic_in_population = sum(allele_type_vec(:)==STOP) / (target_size_by_class_vec(2) * iters);
    empirical_prob_missense_polymorphic_in_population = sum(allele_type_vec(:)==MISSENSE) / (target_size_by_class_vec(3) * iters);    
    
    empirical_prob_neutral_polymorphic_in_sample = sum(sum((X(1:total_num_alleles,:)>0) .* (allele_type_vec==SYNONYMOUS))) / ...
        (target_size_by_class_vec(1) * iters);
    empirical_prob_null_polymorphic_in_sample = sum(sum((X(1:total_num_alleles,:)>0) .* (allele_type_vec==STOP))) / ...
        (target_size_by_class_vec(2) * iters);
    empirical_prob_missense_polymorphic_in_sample = sum(sum((X(1:total_num_alleles,:)>0) .* (allele_type_vec==MISSENSE))) / ...
        (target_size_by_class_vec(3) * iters);

    numeric_prob_neutral_polymorphic_in_sample_given_population = 1-sum(p_neutral_vec .* (1-x_vec).^num_individuals) ./ sum(p_neutral_vec); % Get integral 
    numeric_prob_null_polymorphic_in_sample_given_population = 1-sum(p_null_vec .* (1-x_vec).^num_individuals) ./ sum(p_null_vec); % Get integral 
    
    
    P_poly = var2struct(empirical_prob_neutral_polymorphic_in_population, empirical_prob_neutral_polymorphic_in_sample, ...
        empirical_prob_null_polymorphic_in_population, empirical_prob_null_polymorphic_in_sample, ...
        empirical_prob_missense_polymorphic_in_population, empirical_prob_missense_polymorphic_in_sample, ...
        numeric_prob_neutral_polymorphic_in_sample_given_population, numeric_prob_null_polymorphic_in_sample_given_population); 
else
    P_poly = [];
end
