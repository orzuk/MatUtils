% Determine architecture parameters based on it's type
% 
% Input: 
% N - number of loci 
% architecture_str - type of architecture
% h_interval - constraints on heritability 
% freq_interval - constraints on prevalence 
% penetrance_interval - constraints on penetrance (of what?)
% k_in_clause - number of different SNPs in one pathway 
% num_clauses - number of pathways 
% iters - number of iterations for ...
%
% Output: 
% params_struct - structure with all parameters 
% 
function params_struct = set_architecture_params(N, architecture_str, ...
    h_interval, freq_interval, penetrance_interval, ...
    k_in_clause, num_clauses, iters)


penetrance_interval = [penetrance_interval(1) * freq_interval(2) ...
    penetrance_interval(2) * freq_interval(1)]; % normalize from lambda_mz risk to penetrance
penetrance_interval = [min(penetrance_interval) max(penetrance_interval)]; 

if(~exist('iters', 'var') || isempty(iters))
    iters = 1; 
end
params_struct.linear_coef_vec = rand(N,1); % coefficients of linear combinations
rand_flag = 0; 
if(rand_flag)
    z_std_vec = rand(iters,1); 
else
    z_std_vec = (1:iters) ./ iters;
end
params_struct.z_std = 2 .^ (-5*z_std_vec+ ...
    (1/2).*log2((1-h_interval(1))./4)); % don't waste variances outside h's range. take a different extra variance each time (for binary we can't have std > 0.5)

% Set prevalence (frequency) range. How do we use h information?  Flip min and max!!
params_struct.min_freq_interval = [0 min(freq_interval(2), penetrance_interval(2))]; % Set minimum and maximum values 
params_struct.max_freq_interval = [max(freq_interval(1), penetrance_interval(1)) 1]; % set maximum

params_struct.min_freq = zeros(iters,1); % now loop on iters
params_struct.max_freq = zeros(iters,1); 
grid_size = ceil(sqrt(iters)); % we make a 2D grid 

min_freq_grid = ((0:grid_size-1) ./ (grid_size-1)) .* ...
    (params_struct.min_freq_interval(2)-params_struct.min_freq_interval(1)) + ...
    params_struct.min_freq_interval(1);
max_freq_grid = ((0:grid_size-1) ./ (grid_size-1)) .* ...
    (params_struct.max_freq_interval(2)-params_struct.max_freq_interval(1)) + ...
    params_struct.max_freq_interval(1);
params_struct.min_freq = mat_into_vec(repmat(min_freq_grid, grid_size, 1));
params_struct.max_freq = mat_into_vec(repmat(max_freq_grid', 1, grid_size));
params_struct.max_freq = max(params_struct.min_freq, params_struct.max_freq); % must have a<b
params_struct.min_freq = vec2row(params_struct.min_freq(1:iters));
params_struct.max_freq = vec2row(params_struct.max_freq(1:iters));

if(exist('k_in_clause', 'var'))
    params_struct.k_in_clause = k_in_clause;
end
if(exist('num_clauses', 'var'))
    params_struct.num_clauses = num_clauses;
end

switch architecture_str % prepare the parameters vec
    case { 'sum-of-ands', 'sum-of-ors', 'CNF', 'DNF'}        
        non_overlap_flag = 1; % force clauses to be non-overlapping
        %        k_in_clause = params_vec(N+1); % number of variables appearing in one clause
        %        num_clauses = params_vec(N+2); % number of clauses
        if(non_overlap_flag) % force clauses to be non-overlapping
             params_struct.clause_mat = ...
                 zeros(params_struct.num_clauses, params_struct.k_in_clause);  
            for i=1:params_struct.num_clauses
                params_struct.clause_mat(i,:) = ...
                    (i-1)*params_struct.k_in_clause+1:i*params_struct.k_in_clause;
            end            
        else % choose clauses randomly
            params_struct.clause_mat = nchoosek(1:N, params_struct.k_in_clause); % choose clauses randomly
            C_inds = randperm(size(params_struct.clause_mat,1));
            C_inds = C_inds(1:params_struct.num_clauses);
            params_struct.clause_mat = params_struct.clause_mat(C_inds, :); % this gives indices for the clauses matrix
        end
   
    case {'additive', 'linear'}
        
    case 'random'
        params_struct.linear_coef_vec(1) = 0.1;
    case {'sigmoid', 'sigmoid-additive', 'and-of-sigmoids', 'or-of-sigmoids'} % add a,b,c of sigmoid
%         [params_struct.a params_struct.b params_struct.c] = ...
%             [params_vec' cumsum([1 50 10].*rand(1,3))]';
        params_struct.linear_coef_vec(1:N) = 1; % force all coefficients to be one (simpler model)
        d = sum(params_struct.linear_coef_vec(1:N));
        params_struct.a = rand(1) * d; % set sigmoid in the middle
        params_struct.b = 50*rand(1) / d; % w
        %                    alpha = 0.1*rand(1);
        %                    params_vec{i}(end) = -(1/alpha)*log(params_vec{i}(end-2)/params_vec{i}(end-1)) / sum(params_vec{i}(1:N));
    case 'given'
        params_struct.linear_coef_vec = [0.02 0.01 0 1];
        
    case 'circuit' % New: enable a general circuit 
        
        
    case 'and-of-k_or_more_of_n' % special 
        params_struct.a = 2; % set k - how many are needed in each pathway 
end

params_struct.n_samples_vec = 200:200:10000; % Set power parameters for GWAS 
params_struct.power_alpha = 10^(-8); % level of significance for GWAS
params_struct.power_pairwise_alpha = 10^(-5); % level for epistasis detection should be higher
params_struct.power_iters = 500;
params_struct.power_test_type = 'marginal'; % we don't know how to do epistasis yet
params_struct.power_test_stat = 'chi-square';
params_struct.power_pairwise_test_stat = 'logistic';

