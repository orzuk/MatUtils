% Compute statistics for several architectures
% For now work only for binary architectures.
% Need to change this to include hirerchical functions ('and-of..');
%
% Input:
% architecture_str - type of architecture
% f_vec - minor allele frequencies
% params_struct - various parameters (affine coefficients)
% p_x_vec - Pr(x) for genotype vectors
% p_x_times_z_vec Pr(x, z) for genotype and phenotype (disease)
% iters - number of sampling iterations
% compute_method_flag - how to compute (sample, enumerate or analytic formulas)
%
% Output:
% mu - mean of trait
% V - variance of trait
% v_environment - amount of variance left given all genotypes
% mz_twin_risk - prob. of having trait for an identical twin of a disease person
% H - heritability content
% h_liability - heritability content on liability scale
%
function [mu V v_environment v_genetic mz_twin_risk H h_liability] = ...
    compute_architecture_statistics(architecture_str, ...
    f_vec, params_struct, p_x_vec, p_x_times_z_vec, iters, ...
    compute_method_flag) % compute several moments and other stuff for architecture


AssignGeneralConstants;
AssignStatsConstants;
N = length(f_vec);
q_z = std_to_q_binary(params_struct.z_std);

if( ~exist('x_vec', 'var') || isempty(x_vec)) % initilize x_vec if needed
    switch compute_method_flag
        case 'sampling' % sample a few x's
            %    x_vec = double(rand(iters,N) < 0.5);
            x_vec = double(rand(iters,N) < repmat(f_vec, iters, 1)); % note: used here the same iters for two things !!!
            %        iters = min(100, iters); % lower sampling from multiple x's !!!
        case 'enumerate' % enumerate all 2^N possible vectors
            x_vec =  dec2bin(0:2^N-1) - '0'; % take all 2^N possibilities
        case 'analytic' % analyic computation. No need to enumerate
    end
end
if(~exist('p_x_vec', 'var') || isempty(p_x_vec))
    if(~strcmp(compute_method_flag, 'analytic'))
        normalize_flag = strcmp(compute_method_flag, 'sampling');
        p_x_vec = x_to_prob_vec(x_vec, f_vec, normalize_flag); % Compute probabilities for all x possibilities
    end
end

switch compute_method_flag % how to compute
    case 'enumerate'
        iters = 1; % change meaning of iters
        z = genetic_architecture(x_vec, architecture_str, ...
            params_struct, iters); % generate outputs
        %        V = sum(p_x_vec .* z_clean.^2) / iters - sum(p_x_vec .* z_clean)^2 / iters^2 + params_struct.z_std^2; % variance of trait
        if(~exist('p_x_times_z_vec', 'var') || isempty(p_x_times_z_vec))
            p_x_times_z_vec = vec2row(p_x_vec) .* z; % temporary vector to save time
        end
        mu = sum(p_x_times_z_vec) / iters; % mean of trait (just used for other calculations)
        V = mu*(1-mu);
        v_environment = sum( p_x_times_z_vec .* (1-z)) / iters;
        mz_twin_risk =  (sum(p_x_times_z_vec .* z) / iters) / mu;
        
    case 'sampling'
        
        iters = 1; % change meaning of iters
        z = genetic_architecture(x_vec, architecture_str, ...
            params_struct, iters); % generate outputs
        
        p_x_times_z_vec = vec2row(p_x_vec) .* z; % temporary vector to save time
        mu = sum(p_x_times_z_vec) / iters; % mean of trait (just used for other calculations)
        V = sum(vec2row(p_x_vec) .* vec2row(z)) / iters - mu^2; % variance of trait. z is treated as prob.!!!
        v_environment = sum( p_x_times_z_vec .* (1-z)) / iters;
        mz_twin_risk =  (sum(p_x_times_z_vec .* z) / iters) / mu;
        %        iters = 1; % previously was set here - why???
    case 'analytic' % works for certain (most) architectures
        iters = length(params_struct.z_std); % determine iterations by number of variance parameters
        mu=1; v_environment = 0;
        switch architecture_str
            case 'sigmoid' % we can compute sigmoid analytically only if all the MAF are the same
                %                f = f_vec(1);
                binom_probs = bernoulli_sum_prob(f_vec); %  binopdf(0:N, N, f); % get all binom probs
                mu = sum(binom_probs .* tanh(params_struct.a .* ([0:N] - params_struct.b)));
                v_environment = sum(binom_probs .* ...
                    tanh(params_struct.a .* ([0:N] - params_struct.b)) .* ...
                    (1-tanh(params_struct.a .* ([0:N] - params_struct.b))));
                [mu V v_environment] = ...
                    architecture_affine(mu, v_environment, 0.5, 0.5); % add sigmoid correction for tanh: 0.5 (1 + tanh(a x - b))
            case {'and-of-sigmoids', 'or-of-sigmoids', 'and-of-k-of-n', 'and-of-k_or_more_of_n'} % take a few sigmoids
                switch architecture_str
                    case {'and-of-sigmoids', 'and-of-k-of-n'}
                        mu=1; v_environment = 0;
                    case 'or-of-sigmoids'
                        mu=0; v_environment = 0;
                end
                V=mu*(1-mu);
                pathway_str = architecture_str_to_one_pathway_str(architecture_str);
%                 pathway_delim = strfind(architecture_str, '-');                
%                 architecture_str(pathway_delim(2)+1:end)
                for i=1:params_struct.num_clauses % loop on clauses
                    cur_params_struct = params_struct;
                    cur_params_struct.linear_coef_vec = ...
                        params_struct.linear_coef_vec((i-1)*params_struct.k_in_clause+1:i*params_struct.k_in_clause);
                    cur_params_struct.min_freq(:) = 0; cur_params_struct.max_freq(:) = 1;
                    [mu2 V2 v_environment2] = ...
                        compute_architecture_statistics(pathway_str, ...
                        f_vec((i-1)*params_struct.k_in_clause+1:i*params_struct.k_in_clause), ...
                        cur_params_struct, [], [], iters, ... % we don't need p_x_vec here !!!! 
                        compute_method_flag); % compute several moments and other stuff for architectur
                    switch architecture_str
                        case {'and-of-sigmoids', 'and-of-k-of-n', 'and-of-k_or_more_of_n'}
                            [mu V v_environment] = ...
                                architecture_and(mu, mu2, V, V2, v_environment, v_environment2); % Take AND of two architectures
                        case 'or-of-sigmoids'
                            [mu V v_environment] = ...
                                architecture_or(mu, mu2, V, V2, v_environment, v_environment2); % Take OR of two architectures
                    end
                end
                
                
            case {'CNF', 'DNF', 'sum-of-ors', 'sum-of-ands'}  % we know how to compute for boolean gates
                num_clauses = size(params_struct.clause_mat, 1);
                clause_prob = zeros(num_clauses,1);
                for i=1:num_clauses
                    switch architecture_str
                        case 'CNF' % and of or's
                            clause_prob(i) = 1 - prod(1-f_vec(params_struct.clause_mat(i,:)));
                            mu = mu * clause_prob(i);
                        case 'sum-of-ors'
                            clause_prob(i) = 1 - prod(1-f_vec(params_struct.clause_mat(i,:)));
                            mu = mu + clause_prob(i);
                        case 'DNF' % or of and's
                            clause_prob(i) = 1 - prod(f_vec(params_struct.clause_mat(i,:)));
                            mu = mu * clause_prob(i);
                        case 'sum-of-ands'
                            clause_prob(i) = prod(f_vec(params_struct.clause_mat(i,:)));
                            mu = mu + clause_prob(i);
                    end
                end % loop on number of clauses
                switch architecture_str
                    case 'DNF'
                        mu = 1 - mu;
                    case {'sum-of-ors', 'sum-of-ands'}
                        mu = (mu-1) / num_clauses; clause_prob = clause_prob ./ num_clauses;
                        v_environment =  ...
                            sum(clause_prob .* (1-1/num_clauses)) - ...
                            sum(sum(clause_prob * clause_prob')) + sum(clause_prob.^2); % trick: divided clause_prob instead of p(1-p)
                end
            case {'k-of-n', 'k_of_n', 'k_of_N', 'k_out_of_N'} % this one requires exactly k out of N to be present
                binom_probs = bernoulli_sum_prob(f_vec); %  binopdf(0:N, N, f); % get all binom probs
                mu = binom_probs(params_struct.a + 1); % a represents k (array of probs starts with zero) 
                v_environment = 0; % take a precise architecture (z=1 iff \sum_i x_i = k)
                
            case {'k_or_more_of_n', 'k_or_more_of_N', 'k_out_of_N_or_more', 'at_least_k_of_N'} % this one requires at least k out of N to be present
                binom_probs = bernoulli_sum_prob(f_vec); %  binopdf(0:N, N, f); % get all binom probs
                mu = sum(binom_probs((params_struct.a + 1):end)); % a represents k
                v_environment = 0; % take a precise architecture (z=1 iff \sum_i x_i = k)
                
            case 'circuit' %  a general circuit - circuit propagation not implemented yet .. (go to enumeration)
                error('No analytic architecture statistics for general circuits yet');
        end % switch architecture_str
        if(isfield(params_struct, 'min_freq')) % new: allow any affine transofmation alpha*x+beta (or the other way around ..)
            alpha = params_struct.max_freq - params_struct.min_freq;
            beta = params_struct.min_freq;
        else
            alpha = 1-2*q_z; beta = q_z;
        end
        [mu V v_environment mz_twin_risk] = ...
            architecture_affine(mu, v_environment, alpha, beta); % Take 'AND' of two architectures
        %         v_environment = beta^2 * v_environment + beta*(1-2*alpha-beta) * mu + alpha*(1-alpha); %        v_environment2 = v_environment + params_struct.z_std^2; % mu*(1+mu) - V; % variance unexplained by the genotypes.  % WRONG! (works only for binary clean arch.)
        %         mu = mu*(1-q_z) + (1-mu)*q_z; % convolution with noise
        %         V = mu*(1-mu); % works for any binary variable
        %         mz_twin_risk = 1 - v_environment / mu;
end

v_genetic = V - v_environment; H = v_genetic ./ V; % compute heritability
if( abs(mz_twin_risk ./ mu - ( 1 - H + H ./ mu ) ) > 0.000000001 )
    problem_mz_and_h = 2345325325
    error('mz risk doesnt match heritability');
end




% Simple auxillary functions that help us put together complicated architectures.
% These don't always work like 'regular' boolean gates or random variables
% These functions are used in the analytic computations. They always assume
% that the two (or more) inputs are independent
% Here try to write a GENERAL gate with all options 
function  [mu V v_environment_out mz_twin_risk] = ...
    architecture_gate(mu, v_environment, gate_type, alpha, beta)

mu_out = apply_gate(mu, gate_type, gate_params); % standard: just pass to gate 

switch gate_type
    
    case {'liablility', 'probit'} % compute liability transformation 
        V_environment_out = mu_out - quadl(@(x)broad_fun(x, T, H), -100, 100);         
end
v_environment_out = 0; % ??  


if(nargout > 2) % compute two dependent parameters
    [V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment); % standard
end


% Simple auxillary functions that help us put together complicated architectures.
% These don't always work like 'regular' boolean gates or random variables
% These functions are used in the analytic computations. They always assume
% that the two (or more) inputs are independent
function [mu V v_environment mz_twin_risk] = ...
    architecture_or(mu1, mu2, V1, V2, v_environment1, v_environment2) % Take or of two architectures

[mu V v_environment] = ...
    architecture_and(1-mu1, 1-mu2, V1, V2, v_environment1, v_environment2); % Take and of two architectures
mu = 1-mu; % flip mu 
[V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment); % standard

function [mu V v_environment mz_twin_risk] = ...
    architecture_xor(mu1, mu2, v_environment1, v_environment2) % Take xor of two architectures
mu = mu1.*(1-mu2) + (1-mu1).*mu2;
v_environment = mu - mu1 - mu2 + 5.*v_environment1.*v_environment2 - 14.*mu1.*mu2 + ...
    8.*mu1.*v_environment2 + 8.*mu2.*v_environment1; % Where did this calculation come from ? 
[V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment); % standard

function [mu V v_environment mz_twin_risk] = ...
    architecture_sum(mu1, mu2, v_environment1, v_environment2) % Take sum of two architectures

mu = mu1 + mu2;
v_environment = v_environment1 + v_environment2 - 2.*mu1.*mu2; % wild guess (??)
[V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment);

function [mu V v_environment mz_twin_risk] = ... % works also for a vector of parameters
    architecture_affine(mu1, v_environment1, alpha, beta) % Take affine transformation of an architecture

mu =  alpha.*mu1 + beta;
v_environment = alpha.^2 .* v_environment1 + alpha.*(1-2.*beta-alpha) .* mu1 + beta.*(1-beta);
[V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment);


function [mu V v_environment mz_twin_risk] = ... % works also for a vector of parameters
    architecture_not(mu1, v_environment1) % Take affine transformation of an architecture

[mu V v_environment mz_twin_risk] = architecture_affine(mu1, v_environment1, -1, 0);

function [mu V v_environment mz_twin_risk] = ...
    architecture_k_or_more_of_n(mu1, mu2, v_environment1, v_environment2, k) % Take at least k-of-N architectures

mu = mu1;
%v_environment = v_environment1 + v_environment2 - 2.*mu1.*mu2; % wild guess (??)

[V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment);




