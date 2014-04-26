% Compute many statistics for a disease model 
% with pairwise (and higher-order) interactions
%
% Input:
% f_vec - frequency of each marker
% architecture_str - type of genetic interaction
% params_struct - parameters specific for the genetic architecture % z_std - extra parameter giving the st.d. of noise
% compute_method_flag - says if to sample a few x's or take all 2^N x's (more accurate but very heavy)
% iters - number of instances to generate for each genotype vector
% x_vec - input x genotypes (optional)
% p_x_vec - probabilities of input x genotypes (optional)
%
% Output:
% mu - mean of trait (frequency in population for binary traits, i.e. disease prevalence)
% V - variance of trait
% v_marginal - amount of variance left given each individual locus
% v_environment - amount of variance left given all genotypes
% v_genetic - amount of variance explained by all genotypes
% V_helper - another estimator of V (matching the additive explained estimator)
% v_additive_explained - amount of variance explained by additive model
% v_pairwise - variance of trait given two loci
% v_pairwise_explained - variance explained by two loci
% H - fraction of variance explaied by genetic factors
% h_add - fraction of variance explained by additive genetic factors 
% h_liability - heritability on a liability scale (of fitted model values) 
% mz_twin_risk - prob. of having trait for an identical twin of a someone with disease
% GRR_marginal - ratio of disease prob. between someone with and without risk genotype
% p_z_x_marginal
% lods_ratio_pairwise - extra effect of pairwise interaction on lods ratio (beyond marginal effect)
% mu_given_k_ones - prob. having the disease given k ones in genotype
% mu_given_k_ones_std - prob. having the disease given k ones in genotype - variance
% relative_risk - risk for disease for relatives of different degrees (only direct descendents)
% family_tree - a tree of family members
% family_risk - risk for relatives of a diseased person (all family relations)
% architecture_formula - formula describing the architecture
%
function [candidate_architectures mu V v_marginal v_environment v_genetic V_helper ...
    v_marginal_explained v_additive_explained v_pairwise v_pairwise_explained ...
    H h_add h_liability h_add_fitted h_mult_fitted ...
    H_from_twins h_add_from_twins h_liability_from_twins h_liability_from_twins_ADE ... 
    mz_twin_risk GRR_marginal p_z_x_marginal lods_ratio_pairwise mu_pairwise ...
    mu_given_k_ones mu_given_k_ones_std relative_risk family_tree family_risk ...
    architecture_formula full_calc_flag] = ...
    compute_architecture_statistics_master(f_vec, architecture_str, params_struct, ...
    compute_method_flag, iters, x_vec, p_x_vec, x_ind_mat, ...
    h_interval, h_add_interval, ratio_interval, ...
    freq_interval, lambda_mz_interval, lods_interval)

AssignGeneralConstants;
AssignStatsConstants;
full_calc_flag = 1;
N = length(f_vec); % number of markers
num_x = size(x_vec, 1);
switch compute_method_flag
    case 'sampling'
        p_x_vec(:) = 1 / num_x;
        iters = min(iters, 100);
    case 'analytic'
        iters = length(params_struct.z_std); % enable vectorization
    otherwise
        iters = 1; % for all x's we need to randomize only once
end

v_marginal = zeros(N,iters); GRR_marginal = zeros(N,iters);
p_z_x_marginal = zeros(N,4,iters); % Prob (z=i,x=j)
v_pairwise = zeros(N,N,iters); lods_ratio_pairwise = zeros(N,N,iters);  % compute pairwise interactions
mu_pairwise = zeros(N,N,iters,4); % Pr(z=1|x_i,x_j)
v_pairwise_explained = zeros(N,N,iters); % Var(z) - Var(z|x_i,x_j) ??? 
v_marginal_explained = zeros(N,iters); % Var(z) - Var(z|x_i)
v_additive_explained = sum(v_marginal_explained); % N*Var(z) - \sum_j Var(z|x_j)
h_liability = zeros(1,iters); % liability 
h_add_fitted = zeros(1,iters); % additive fitted
h_mult_fitted = zeros(1,iters); % multiplicative fitted 
H_from_twins = zeros(1,iters); % estimated broad-sense from MZ twins relative-risk
h_add_from_twins = zeros(1,iters); % estimated narrow-sense from MZ and DZ twins
h_liability_from_twins = zeros(1,iters); % estimated  liability-threshold from MZ and DZ twins (ACE model) 
h_liability_from_twins_ADE = zeros(1,iters); % estimated  liability-threshold from MZ and DZ twins (ADE model)

ttt_arch = cputime;

switch compute_method_flag
    case {'sampling', 'enumerate'}
        z  = genetic_architecture(x_vec, architecture_str, ...
            params_struct, iters); % generate outputs
    otherwise
        z = 1;
end
% ttt_arch = cputime - ttt_arch;
% ttt_start = cputime;

arch_type = architecture_to_type(architecture_str);
architecture_formula = get_architecture_formula(architecture_str, params_struct, 1);
switch compute_method_flag
    case 'sampling'
        mu_helper = sum(p_x_vec .* z(1:num_x));
        V_helper =  sum(p_x_vec .* z(1:num_x).^2) - mu_helper^2; % compute based on only first set
    otherwise
        V_helper = [];
end
p_x_vec = repmat(p_x_vec, iters, 1); % duplicate to match iters
p_x_times_z_vec = p_x_vec .* z; % temporary vector to save time

% compute_marginals_flag = 1; 
compute_pairwise_flag = 1;
[mu V v_environment v_genetic mz_twin_risk H] = ...
    compute_architecture_statistics(architecture_str, ...
    f_vec, params_struct, p_x_vec, p_x_times_z_vec, iters, ...
    compute_method_flag); % compute several moments and other stuff for architecture
lambda_mz = mz_twin_risk ./ mu; 
full_calc_flag = (mu >= freq_interval(1)) & (mu <= freq_interval(2));
full_calc_flag = full_calc_flag & ((lambda_mz >= lambda_mz_interval(1)) & ...
    (lambda_mz <= lambda_mz_interval(2)));
full_calc_flag = full_calc_flag & ((H >= h_interval(1)) & (H <= h_interval(2)));
full_calc_inds = find(full_calc_flag);
if(any(H(full_calc_inds) < -sqrt(epsilon)) || any(H(full_calc_inds) > 1+epsilon)) % something's wrong
    ttt_h = 9999
end

switch compute_method_flag
    case 'sampling'
        iters = 1;
end
if(~isempty(full_calc_inds))
    save_z_std = params_struct.z_std;
    params_struct.z_std = params_struct.z_std(full_calc_inds); % take only relevant architectures
    if(isfield(params_struct, 'min_freq'))
        save_min_freq = params_struct.min_freq;
        save_max_freq = params_struct.max_freq;
        params_struct.min_freq = params_struct.min_freq(full_calc_inds);
        params_struct.max_freq = params_struct.max_freq(full_calc_inds);
    end
    [v_marginal(:,full_calc_inds) GRR_marginal(:,full_calc_inds) ...
        p_z_x_marginal(:,:,full_calc_inds)] = ...
        compute_architecture_statistics_marginal(architecture_str, ...
        f_vec, params_struct, p_x_vec, p_x_times_z_vec, iters, ...
        compute_method_flag, mu(full_calc_inds)); % compute several moments and other stuff for architecture
    v_marginal_explained = repmat(V, N, 1) - v_marginal;
    v_additive_explained = sum(v_marginal_explained);
    h_add = v_additive_explained ./ V; r = h_add ./ H;
    
    
    full_calc_flag = full_calc_flag & ((h_add >= h_add_interval(1)) & ...
        (h_add <= h_add_interval(2)));
    full_calc_flag = full_calc_flag & ((r >= ratio_interval(1)) & ...
        (r <= ratio_interval(2)));
    full_calc_flag = full_calc_flag & ((max(GRR_marginal) >= lods_interval(1)) & ...
        (max(GRR_marginal) <= lods_interval(2)));
end

%ttt_pairwise = cputime;
full_calc_inds = find(full_calc_flag);
if(~isempty(full_calc_inds)) % full_calc_flag)
    params_struct.z_std = save_z_std (full_calc_inds); % take only relevant architectures
    if(isfield(params_struct, 'min_freq'))
        params_struct.min_freq = save_min_freq(full_calc_inds);
        params_struct.max_freq = save_max_freq(full_calc_inds);
    end
    
    [v_pairwise(:,:,full_calc_inds) v_pairwise_explained(:,:,full_calc_inds) ...
        lods_ratio_pairwise(:,:,full_calc_inds) mu_pairwise(:,:,full_calc_inds,:)] = ...
        compute_architecture_statistics_pairwise(architecture_str, ...
        V(full_calc_inds), f_vec, params_struct, p_x_vec, p_x_times_z_vec, x_ind_mat, iters, ...
        v_marginal(:,full_calc_inds), GRR_marginal(:,full_calc_inds), ...
        compute_method_flag, compute_pairwise_flag); % compute several moments and other stuff for architecture
end



sampling_iters = 100; % just for computing mu given k
%mu_given_k_ones = zeros(N+1,iters); mu_given_k_ones_std = zeros(N+1,iters);
mu_given_k_ones = zeros(N,N+1,iters); mu_given_k_ones_std = zeros(N,N+1,iters); % new: get matrices 
mu_given_k_one_pathway = zeros(params_struct.k_in_clause, params_struct.k_in_clause+1, iters); 
max_generations = 1; 
relative_risk = zeros(max_generations, iters); % saves risk for a few generations
family_tree = generate_family_tree(max_generations); family_size = length(family_tree); 
family_risk = zeros(family_size, iters); 
if(~isempty(full_calc_inds))
    params_struct.z_std = save_z_std (full_calc_inds);
    if(isfield(params_struct, 'min_freq'))
        params_struct.min_freq = save_min_freq(full_calc_inds);
        params_struct.max_freq = save_max_freq(full_calc_inds);
    end
    
    %        mu2(full_calc_inds) v_environment2(full_calc_inds)] = ...
    [mu_given_k_ones(:,:,full_calc_inds) mu_given_k_ones_std(:,:,full_calc_inds) ...
        mu_given_k_one_pathway(:,:,full_calc_inds) ...
        relative_risk(:,full_calc_inds) family_tree family_risk(:, full_calc_inds) ... 
        h_liability(full_calc_inds) h_add_fitted(full_calc_inds) h_mult_fitted(full_calc_inds) ...
                H_from_twins(full_calc_inds) h_add_from_twins(full_calc_inds) ...
                h_liability_from_twins(full_calc_inds) h_liability_from_twins_ADE(full_calc_inds)] = ...
        compute_architecture_statistics_extra(architecture_str, ...
        f_vec, params_struct, sampling_iters, ...
        mu(full_calc_inds), V(full_calc_inds), h_add(full_calc_inds), ...
        v_environment(full_calc_inds), GRR_marginal(:,full_calc_inds), ... % add architecture parameters as input
        compute_method_flag, max_generations); % compute other stuff for architecture (this is not very time-costly)
end
H = v_genetic ./ V; h_add = v_additive_explained ./ V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
candidate_architectures = []; % finished entire function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isnan(v_additive_explained )) % Testing of unexpected results
    ttt = 999
end
switch compute_method_flag
    case 'sampling'
    otherwise
        if(max(v_marginal) > V + epsilon) % marginal larger than
            ttt = 999;
        end
        if( (~isempty(v_additive_explained > V) || ~isempty(v_additive_explained < 0))) % too much explained ... not possible
            ttt = 999;
        end
end
% if(~isempty(v_additive_explained ./ V < 0.2)) % wrong only for a specific model
%     ttt = 999;
% end
% if(~isempty(v_additive_explained ./ v_genetic < 0.2)) % wrong only for a specific model
%     ttt = 999;
% end
if(any(v_additive_explained(full_calc_inds) > v_genetic(full_calc_inds) + sqrt(epsilon)))
    error_is = 'Error! additive var > total explained var! cant be!!!'
    ttt = 999988
end
if(any(v_genetic(full_calc_inds) > V(full_calc_inds) + sqrt(epsilon)) || ...
        any(v_genetic(full_calc_inds) < -sqrt(epsilon)))
    tttt = 999999
end
if(any(v_additive_explained(full_calc_inds) < -sqrt(epsilon)))
    ttt = 124
end
