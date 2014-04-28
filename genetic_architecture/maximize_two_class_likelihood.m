% Compute the maximum likelihood estimator given genotype and phenotype data for the two class model
%
% Input:
% s_null_vec - vector of possible selection coefficients for the null alleles. Should be NEGATIVE for delterious alleles (so fitness is 1+s)
% alpha_vec - vector of possible fraction of null alleles at birth
% beta_vec - vector of possible effect size of null alleles
% rare_cumulative_per_gene - total expected number of rare alleles in gene (related to theta)
% target_size_by_class_vec - NEW! target size for each allele type (used in poisson model). This is like the above but counts expected number of different alleles, not their frequencies !  %%% NEW! alleles class type vector. This states for each allele if it is neutral, harmfull, or in a mixture (missense) 
% N - effective population size
% X - genotype data matrix
% y - phenotype data vector
% trait_type - disease or quantitative
% prevalence - frequency of disease for disease trait
% null_w_vec - (optional) this is the assignment of which alleles are null
%            and which not. This simplifies the likelihood computation alot when it is known
% maximize_parameters - on which parameters to do the maximization. Default is all (s, alpha, beta). Format is [S alpha beta]  (all binary variables)
% full_flag - how is input data (X and y) given 
% D - demographic model which controls the allele frequncy distribution (NEW! Default is constant population size)
% num_individuals - 
% implementation_str how to find minimum 
% 
% Output:
% max_LL - maximum of log-likelihood of data
% max_alpha - value of alpha attaining maximum
% max_beta - value of beta attaining maximum
% max_s - value of s attaining maximum
% NEW! should add also confidence intervals (not implemented yet) (Use Fisher's information matrix)  
%
function [max_LL max_s max_alpha max_beta] = ...
    maximize_two_class_likelihood(s_null_vec, alpha_vec, beta_vec, rare_cumulative_per_gene, target_size_by_class_vec, N, ...
    X, y, trait_type, prevalence, null_w_vec, maximize_parameters, full_flag, D, num_individuals, implementation_str)


if(~exist('implementation_str', 'var') || isempty(implementation_str))
    implementation_str = 'minsearch'; % 'brute-force'; % how to find optimum 
end
    
if(~exist('trait_type', 'var') || isempty(trait_type))
    trait_type = [];
end
if(~exist('prevalence', 'var') || isempty(prevalence))
    prevalence = [];
end
if(~exist('null_w_vec', 'var') || isempty(null_w_vec))
    null_w_vec = [];
end
if(~exist('maximize_parameters', 'var') || isempty(maximize_parameters))
    maximize_parameters = [1 1 1]; % maximize over all parameters
end
if(~exist('target_size_by_class_vec', 'var') || isempty(target_size_by_class_vec))   
    poisson_model_flag = 0;  expand_format_flag = 'individual'; 
else
    poisson_model_flag = 1; expand_format_flag = 'summary'; % use poisson model for each class 
end
if(~exist('include_phenotype', 'var') || isempty(include_phenotype)) % default: don't include phenotypes 
    % include_phenotype = 0; 
    include_phenotype = maximize_parameters(3); % include phenotypes if and only if we optimze over beta 
end

if(~exist('full_flag', 'var') || isempty(full_flag))
    full_flag = 1; 
end
if(full_flag)
    [num_individuals L] = size(X); % set number of individuals and number of SNPs
else
    if(~isempty(y)) % when phenotype is given
        num_individuals = length(y);
        [~, ~, ~, ~, L] = ... % get numbers from compact representation
            expand_two_class_summary_statistics(X, num_individuals, expand_format_flag);
    else % here we have no way of knowing the # of individuals !!! assume they're given in advance (BUT! might be different number for different alleles!)
        L = length(null_w_vec); % get number of alleles 
    end
end

switch implementation_str  % choose how to maximize likelihood
    case 'brute-force'  % here simply enumerate grid-points and find point with ML
        if(~maximize_parameters(3)) % don't loop over beta
            beta_vec = 0; % meaningless
        end
        % PROBLEM: HERE LOG-LIKELIHOOD IS MONOTONICALLY DECREASING WITH
        % ALPHA - WHY? 
        log_like_mat = ... % compute likelihood (currently vary only alpha)
            compute_two_class_log_likelihood(s_null_vec, alpha_vec, beta_vec, rare_cumulative_per_gene, target_size_by_class_vec, N, ...
            X, y, trait_type, prevalence, null_w_vec, include_phenotype, full_flag, num_individuals);        
        [max_LL tmp_ind] = max(log_like_mat(:)); % find the maximal grid-point
        [I J K] = ind2sub(size(log_like_mat), tmp_ind); 
        max_s = s_null_vec(I); 
        max_alpha = alpha_vec(J); 
        max_beta = beta_vec(K); 
        
    case 'gradient' % use gradient-descent algorithm (??)
        % use gradient descent. First optimize only s and alpha using
        % genotypes. Then optimize beta using phenotype
        iters = 1000; % number of gradient iterations
        
        
    case 'minsearch' % use Matlab's optimization
        if(isempty(null_w_vec) || poisson_model_flag) % here we don't know alpha and the null positions
            [max_s_alpha max_LL_genotype] = fminsearch(@(s_alpha) ... % use genotypes
                -compute_two_class_log_likelihood(s_alpha(1), s_alpha(2), [], rare_cumulative_per_gene, target_size_by_class_vec, ...
                N, X, [], [], [], null_w_vec, include_phenotype, full_flag, num_individuals, D), ...
                [s_null_vec; alpha_vec]); % here s_null_vec and alpha_vec serve as initial guesses! (should be scalars)
            % use phenotypes            
            max_s = max_s_alpha(1); max_alpha = max_s_alpha(2); % We need to check that 0 < alpha < 1
            
        else % here we do know alpha
            [max_s max_LL_genotype] = fminsearch(@(s_alpha) ... % use genotypes
                -compute_two_class_log_likelihood(s_alpha, 1, [], rare_cumulative_per_gene, target_size_by_class_vec, N, X, [], ...
                trait_type, prevalence, null_w_vec, include_phenotype, full_flag, num_individuals), s_null_vec); % maximize over genotypes only 
            max_alpha = 1;
                        
            if(maximize_parameters(3)) % optimize also on effect size beta - here we must use phenotype !! (slower)
                [max_beta max_LL_phenotype] = fminsearch(@(beta) ... % now optimize only phenotype part
                    -compute_two_class_log_likelihood(max_s, 1, beta, rare_cumulative_per_gene, target_size_by_class_vec, N, X, y, ...
                    trait_type, prevalence, null_w_vec, -1, full_flag), beta_vec);
            else
                max_LL_phenotype = ... % Compute likelihood for beta=0
                    compute_two_class_log_likelihood(max_s, 1, 0, rare_cumulative_per_gene, target_size_by_class_vec, N, X, y, ...
                    trait_type, prevalence, null_w_vec, -1, full_flag);
                max_beta = 0;
            end
            max_LL = max_LL_genotype + max_LL_phenotype;            
        end % if null_w_vec is empty (i.e. we don't know alpha)         
end % switch method of MLE estimation


