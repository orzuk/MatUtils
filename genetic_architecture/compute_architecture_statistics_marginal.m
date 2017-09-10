% Compute marginal statistics for several architectures
% This works for now only for binary architectures
%
% Input:
% architecture_str - string representing architecture type
% f_vec - frequencies of each locus
% params_struct - general architecture parameters
% p_x_vec - probability of each genotype vector
% p_x_times_z_vec - probability of each genotype vector times prob. phenotype
% iters - number of iterations for sampling (? )
% compute_method_flag - which method to use (sampling/enumerate/analytic)
% mu - trait's mean (prevalence) (seems to be not needed).
%
% Output:
% v_marginal - vector of variances left in disease given each marker, Var(Z | x_i)
% GRR_marginal - genetic relative risk for each marker: Pr(Z=1|x_i=1) / Pr(Z=1|x_i=0)
% p_z_x_marginal - joint probability of each x_i and z. A 2x2 table represented as a vector of length 4 for each locus
% h_add - additive heritability
%
function [v_marginal GRR_marginal p_z_x_marginal h_add h_liability] = ...
    compute_architecture_statistics_marginal(architecture_str, ...
    f_vec, params_struct, p_x_vec, p_x_times_z_vec, iters, ...
    compute_method_flag, mu) % compute several moments and other stuff for architecture

AssignGeneralConstants;
AssignStatsConstants;
N = length(f_vec); % number of loci
v_marginal = zeros(N,1); GRR_marginal = zeros(N,1);

if( ~exist('x_vec', 'var') || isempty(x_vec)) % initilize x_vec if needed (always needed ..)
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
    normalize_flag = strcmp(compute_method_flag, 'sampling');
    switch compute_method_flag
        case 'sampling'
            p_x_vec = x_to_prob_vec(x_vec, repmat(0.5, 1, N), normalize_flag); % Compute probabilities for all x possibilities
        case 'enumerate'
            p_x_vec = x_to_prob_vec(x_vec, f_vec, normalize_flag); % Compute probabilities for all x possibilities
    end
end


arch_type = architecture_to_type(architecture_str);
switch compute_method_flag % how to compute
    case {'enumerate', 'sampling'}
        iters=1; % only one z for each input vector  
        z = genetic_architecture(x_vec, architecture_str, ...
            params_struct, iters); % generate outputs
        p_z_x_marginal = zeros(N,4);
        
        %        if(compute_marginals_flag)
        for i=1:N % compute variance given each one of the genotypes seperately
            ind_one_vec = find(x_vec(:,i)); num_ones = length(ind_one_vec);
            ind_one_vec = repmat(ind_one_vec, iters, 1); num_x = size(x_vec, 1);
            for j=2:iters % correct indices according to iterations
                ind_one_vec(((j-1)*num_ones+1):(j*num_ones)) = ...
                    ind_one_vec(((j-1)*num_ones+1):(j*num_ones)) + num_x * (j-1);
            end
            if(~exist('p_x_times_z_vec', 'var') || isempty(p_x_times_z_vec))
                p_x_times_z_vec = vec2row(p_x_vec) .* z; % temporary vector to save time
            end
            switch compute_method_flag
                case 'sampling'
                    ind_zero_vec = setdiff([1:num_x]', ind_one_vec);
                case 'enumerate'
                    ind_zero_vec = setdiff([1:length(p_x_vec)]', ind_one_vec);
            end
            f_vec_empirical = sum(p_x_vec(ind_one_vec)) / iters; % Prob(x_i = 1);
            switch arch_type
                case CONTINUOUS % also applies to probabilities (?)
                    if(f_vec_empirical > 0)
                        v_marginal(i) = v_marginal(i) +  (f_vec(i) / f_vec_empirical) * ...
                            ( sum( (p_x_vec(ind_one_vec) ./ iters).* (z(ind_one_vec) .^ 2) ) - ...
                            sum((p_x_vec(ind_one_vec) ./ iters).* z(ind_one_vec)) .^ 2 ./ f_vec_empirical ); %%% + ... % ((1-f_vec(i)) / (1-f_vec_empirical)) *
                    end
                    if(f_vec_empirical < 1)
                        v_marginal(i) = v_marginal(i) + ((1-f_vec(i)) / (1-f_vec_empirical)) * ...
                            ( sum( (p_x_vec(ind_zero_vec) ./ iters) .* z(ind_zero_vec) .^ 2 ) - ...
                            sum((p_x_vec(ind_zero_vec) ./ iters).* z(ind_zero_vec)) .^ 2 ./ (1-f_vec_empirical) );
                    end
                case BINARY
                    tmp_prob_one = sum( (p_x_times_z_vec(ind_one_vec) ./ iters) ) ...
                        / f_vec_empirical; % Pr(z=1 | x_i=1)
                    tmp_prob_zero = sum( (p_x_times_z_vec(ind_zero_vec) ./ iters) ) ...
                        / (1-f_vec_empirical); % Pr(z=1 | x_i=0)
                    v_marginal(i) = f_vec_empirical * tmp_prob_one * (1-tmp_prob_one) + ...
                        (1-f_vec_empirical) * tmp_prob_zero * (1-tmp_prob_zero);
                    GRR_marginal(i) = tmp_prob_one / tmp_prob_zero;
                    
                    p_z_x_marginal(i,1) = (1-f_vec_empirical) * (1-tmp_prob_zero); % Pr(x_i=0,z=0)
                    p_z_x_marginal(i,2) = (1-f_vec_empirical) * tmp_prob_zero; % Pr(x_i=0,z=1)
                    p_z_x_marginal(i,3) = f_vec_empirical * (1-tmp_prob_one); % Pr(x_i=1,z=0)
                    p_z_x_marginal(i,4) = f_vec_empirical * tmp_prob_one; % Pr(x_i=1,z=1)
            end
            ind_vec = zeros(N,1); ind_vec(i) = 1;
            %    v_marginal2(i) = conditional_var(x_vec, p_x_vec, p_x_times_z_vec, ind_vec, iters);
        end % loop on marginals
        %       end % if compute marginals
        
        switch arch_type
            case CONTINUOUS
                switch compute_method_flag
                    case 'sampling'
                        
                    otherwise
                        v_marginal = v_marginal + params_struct.z_std^2;
                end
        end
        
        %     case 'sampling' % do nothing?
        %         iters = 1; p_z_x_marginal = zeros(N,4); % first set iters to one
        %         z = genetic_architecture(x_vec, architecture_str, ...
        %             params_struct, iters); % generate outputs
        
        
    case 'analytic' % works for certain architectures
        iters = length(params_struct.z_std);
        %        if(compute_marginals_flag) % depends on recursion level
        mu_marginal = cell(2,1); v_marginal = cell(2,1); v_environment_marginal = cell(2,1);
        for j=0:1 % why two types of marginals? 
            mu_marginal{j+1} = zeros(N,iters);
            v_marginal{j+1} = zeros(N,iters);
            v_environment_marginal{j+1} = zeros(N,iters);
        end
        %            mu_marginal = zeros(N,2); v_marginal = zeros(N,2);
        %            v_environment_marginal = zeros(N,2);
        for i=1:N % loop on individual variables and use recursion to compute their individual contributions to variance 
            f_vec_marginal = f_vec;
            for j=0:1
                f_vec_marginal(i) = j; % set one variable to either zero or one 
                [mu_marginal{j+1}(i,:) v_marginal{j+1}(i,:) v_environment_marginal{j+1}(i,:)] = ...
                    compute_architecture_statistics(architecture_str, ...
                    f_vec_marginal, params_struct, p_x_vec, p_x_times_z_vec, iters, ...
                    compute_method_flag); % compute several moments and other stuff for architecture
            end
        end
        GRR_marginal = mu_marginal{2} ./ mu_marginal{1}; % Pr(z=1|x=1) / Pr(z=1|x=0)
        v_marginal = v_marginal{1} .* (1-repmat(vec2column(f_vec), 1, iters)) + ...
            v_marginal{2} .* repmat(vec2column(f_vec), 1, iters); % average over values of variable 
        p_z_x_marginal = zeros(N,4,iters); % Prob (z=i,x=j)
        p_z_x_marginal(:,1,:) = repmat((1-vec2column(f_vec)), 1, iters) .* (1-mu_marginal{1});
        p_z_x_marginal(:,2,:) = repmat((1-vec2column(f_vec)), 1, iters) .* mu_marginal{1};
        p_z_x_marginal(:,3,:) = repmat(vec2column(f_vec), 1, iters) .* (1-mu_marginal{2});
        p_z_x_marginal(:,4,:) = repmat(vec2column(f_vec), 1, iters) .* mu_marginal{2};
        %        end
end

%if(~exist('mu', 'var') || isempty(mu)) % compute mu ourselves (should be better!)
if(exist('p_x_times_z_vec', 'var') && ~isempty(p_x_times_z_vec))
    mu = sum(p_x_times_z_vec);
end
V = mu.*(1-mu); % fast&dirty computation of mu (not needed)
bad_marginals = find(v_marginal > repmat(V, N, 1)); % can't have marginal variance larger than total variance since we condition on x_i
if(~isempty(bad_marginals))
    xxx = 999; 
end
v_marginal_explained = repmat(V, N, 1) - v_marginal;
v_additive_explained = sum(v_marginal_explained);
h_add = v_additive_explained ./ V;

if(nargout > 4) % output also liability
    h_liability = zeros(size(h_add)); 
    for i = 1:length(h_add)
        h_liability(i) = heritability_scale_change(h_add(i), 'liability', mu(i)); % Different: just use scale-change
    end
end
