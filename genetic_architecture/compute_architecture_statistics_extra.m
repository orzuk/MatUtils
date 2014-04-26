% Compute several extra statistics for several architectures.
% For now work only for binary architectures
%
% Input:
% architecture_str - string representing the type of architecture
% f_vec - MAF of different loci
% params_struct - specific parameters of architecture
% sampling_iters - number of iterations for computing ... via simulations
% mu - arch. mean
% V - arch. variance
% v_environment - arch. environmental variance
% compute_method_flag - how to compute (analytic, sampling or full enumeration)
% max_generations - for how many generations to compute the risk
%
% Output:
% mu_given_k_ones - mean of trait prob. given k ones
% mu_given_k_ones_std - std of trait prob. given k ones
% mu_given_k_one_pathway - the same but for a given pathway
% relative_risk - prob. having the disease given that a relative has it
% familty_tree - a tree of family members
% family_risk - risk for various family members Pr(z=1|z_relative = 1) (not just direct descendents)
% h_liability - heritability on an unobsered liability scale (the 'true' value)
% h_from_twins - heritability computed from twins
% h_liability_from_twins - heritability computed form twins on the liability scale
%
function [mu_given_k_ones mu_given_k_ones_std mu_given_k_one_pathway ...
    relative_risk family_tree family_risk ...
    h_liability h_add_fitted h_mult_fitted ...
    H_from_twins h_add_from_twins h_liability_from_twins h_liability_from_twins_ADE] = ...
    compute_architecture_statistics_extra(architecture_str, ...
    f_vec, params_struct, sampling_iters, ...
    mu, V, h_add, v_environment, GRR_marginal, ... % add architecture parameters as input
    compute_method_flag, max_generations) % compute several moments and other stuff for architecture

AssignGeneralConstants;
AssignStatsConstants;
N = length(f_vec);
if(~exist('max_generations', 'var') || isempty(max_generations))
    max_generations = 3; % how many generations should we compute
end
relative_sampling_iters = 50*sampling_iters; % used to be 500. Went down to 50 as it's way slower now (so less accurate) ...
% [family_risk family_tree] = compute_architecture_family_risk(); % have a seperate function for family risk

switch compute_method_flag % how to compute
    case 'enumerate'
    case {'analytic', 'sampling'} % for now we know how to copmute relative risk only using sampling
        
        %    case 'analytic' % works for certain architectures
        %        sampling_iters = iters; % how many sampling of permuations are needed
        iters = length(params_struct.z_std); % determine iterations by number of variance parameters
        mu_given_k_ones = zeros(N,N+1,iters); mu_given_k_ones_std = zeros(N,N+1,iters);
        % m_given_k_ones (i,j,k) is Prob.(z = 1 | j-1 x's of the first i are one) in iteration #k
        
        
        rand_perms = zeros(relative_sampling_iters, N);
        for i=1:relative_sampling_iters
            rand_perms(i,:) = randperm(N);
        end
        for M=N:N % loop of number of variables observed (entire loop is quadratic :-( )
            %            cur_f_vec = f_vec; % just copy frequencies
            for i=1:M+1 % loop on number of variables set to one out of M (so we observe i-1 ones, from 0 to M)
                mu_given_k_ones_tmp = zeros(relative_sampling_iters, iters);
                cur_f_vec = f_vec; cur_f_vec(1:M) = zeros(M,1); cur_f_vec(1:i-1) = 1; % set i-1 of the first M to one
                %                cur_f_vec = cur_f_vec(rand_perms); % randomize where the positions of ones are
                %             for j=1:sampling_iters % loop on different samples
                %                 mu_given_k_ones_tmp(j,:) = ...
                %                     compute_architecture_statistics(architecture_str, ...
                %                     cur_f_vec(j,:), params_struct, p_x_vec, p_x_times_z_vec, iters, ...
                %                     compute_method_flag); % compute several moments and other stuff for architecture
                %             end
                cur_x_vec = rand(relative_sampling_iters,N) < cur_f_vec(rand_perms);
                for j=1:iters % loop on different architectures
                    cur_params_struct = params_struct; cur_params_struct.z_std = params_struct.z_std(j);
                    if(isfield(cur_params_struct, 'min_freq'))
                        cur_params_struct.min_freq = params_struct.min_freq(j);
                        cur_params_struct.max_freq = params_struct.max_freq(j);
                    end
                    mu_given_k_ones_tmp(:,j) = ... % Different: just compute the architeture values (much faster)
                        genetic_architecture(cur_x_vec, architecture_str, ...
                        cur_params_struct, 1);
                end
                mu_given_k_ones(M,i,:) = mean(mu_given_k_ones_tmp);
                mu_given_k_ones_std(M,i,:) = std(mu_given_k_ones_tmp);
            end % loop on i, number of loci set to one
        end % loop on M, number of loci observed
        
        
        for M=1:N-1 % Alternative: we've computed for one (M=N), so now just need to compute it for the partial observations
            cur_binom_vec = bernoulli_sum_prob(f_vec(M+1:N)); % (alternative works for mean mu but currently we don't know how to make it work for std.)
            for i=1:M+1 % loop on number of variables set to one out of M (so we observe i-1 ones, from 0 to M)
                mu_given_k_ones(M,i,:) = ...
                    sum( repmat(cur_binom_vec, iters, 1)' .* shiftdim(mu_given_k_ones(N,i:i+N-M,:), 1) );
            end
        end
        %         for K = 1:N
        %             p_vec = binopdf(0:K, K, 0.2); mumumu(K) = sum(p_vec .* mu_given_k_ones(K,1:K+1,4));
        %         end
        %         figure; hold on; plot(1:N, mumumu, '.'); title('estimators for \mu');
        %         plot(0:N, repmat(mu(4), N+1,1), 'r');
        %         K=12; figure; plot(0:K, mu_given_k_ones(K,1:K+1,4))
        
        %         for i=1:N+1 % alternative: loop on number of variables set to one out of N (so we observe i-1 ones, from 0 to M)
        %             params_struct_k_of_N = params_struct;
        %             params_struct_k_of_N.arch = 'k_of_N'; params_struct_k_of_N.a = i-1;
        %              [mu_k_of_N V_k_of_N v_environment_k_of_N v_genetic_k_of_N penetrance_k_of_N] = ...
        %                 compute_architecture_statistics('k_of_N', ...
        %                 f_vec, cur_params_struct, [], [], iters, ...
        %                     'analytic');
        %
        %                 [mu_given_k_ones(M,i,:) mu_given_k_ones_std(M,i,:) ...
        %                     v_environment_give_k_ones penetrance_given_k_ones] = ...
        %                     architecture_and(mu, mu_k_of_N, V, V_k_of_N, v_environment, v_environment_k_of_N) % Take and of two architectures
        %         end
end % switch compute methods

cur_params_struct = params_struct; % take first clause
cur_params_struct.linear_coef_vec = ...
    params_struct.linear_coef_vec(1:params_struct.k_in_clause);
one_pathway_str = architecture_str_to_one_pathway_str(architecture_str);
if(isempty(strmatch(one_pathway_str, architecture_str, 'exact')))
    mu_given_k_one_pathway = ...
        compute_architecture_statistics_extra(one_pathway_str, ...
        f_vec(1:params_struct.k_in_clause), cur_params_struct, sampling_iters, ...
        mu, V, h_add, v_environment, GRR_marginal, ... % add architecture parameters as input
        compute_method_flag, max_generations);
end


if(nargout > 3) % compute only as needed
    [family_risk family_tree relative_risk sibs_genotype sibs_freqs ...
        H_from_twins h_add_from_twins h_liability_from_twins h_liability_from_twins_ADE ...
        lambda_mz_from_twins lambda_dz_from_twins] = ...
        compute_architecture_family_risk(architecture_str, ...
        f_vec, params_struct, relative_sampling_iters, ...
        compute_method_flag, max_generations);
    
    fit_mode = 'liability'; % New: fit a different architecture to match this architecture
    num_architectures = length(mu);
    for i=1:num_architectures
        cur_params_struct = params_struct;
        cur_params_struct.min_freq = params_struct.min_freq(i);
        cur_params_struct.max_freq = params_struct.max_freq(i);
        cur_params_struct.z_std = params_struct.z_std(i);
        %         if(isfield(params_struct, 'clause_mat'))
        %             cur_params_struct.clause_mat = cur_params_struct.clause_mat(:,i);
        %         end
        
        [alpha beta fitted_params_struct] = ...
            fit_architecture_parameters(cur_params_struct, ...
            architecture_str, f_vec, fit_mode, [], compute_method_flag); % GRR marginal not used here
        %        h_liability(i) = sum(alpha.^2 .* vec2column(f_vec.*(1-f_vec))); % fraction of variance explained by genetics on liability scale
        h_liability(i) = heritability_scale_change(h_add(i), 'liability', mu(i)); % Different: just use scale-change
        h_add_fitted(i) = sum(alpha); % heritability on binary scale
        
        
        [lambda_s_vec_fitted ...
            lambda_s_add_fitted lambda_mz_add_fitted h_add_fitted V_add_fitted ...
            lambda_s_mult_fitted lambda_mz_mult_fitted h_mult_fitted V_mult_fitted ...
            sib_freq_mat] = ...
            genetic_relative_risk_to_heritability(f_vec, GRR_marginal(:,i), ...
            mu(i));
    end
    
end
