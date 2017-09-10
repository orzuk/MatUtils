% Compute risk for family members of a diseased person \lambda_R
% Currently method is via sampling
%
% Input:
% architecture_str - name of architecture
% f_vec - single locus frequencies
% params_struct - global architecture parameters
% relative_sampling_iters - number of iterations for sampling relatives
% compute_method_flag - how to compute (sampling/analytic)
% max_generations - size of family tree
%
% Output:
% family_risk - vector of risk for each family member
% family_tree - tree of family members
% relative_risk - again vector of risks (probably not needed)
% sibs_genotype - genotype vector for sibs
% sibs_freq - joint frequency of disease in person and sibs ???
% H_from_twins - estimate of (broad-sense) heritability from twins statistics
% h_add_from_twins - estimate of (narrow-sense) heritability from twins statistics%
% h_liability_from_twins - estimate of (liability-threshold) heritability from twins statistics%
%
function [family_risk family_tree relative_risk sibs_genotype sibs_freqs ...
    H_from_twins h_add_from_twins h_liability_from_twins h_liability_from_twins_ADE ...
    lambda_mz lambda_dz] = ...
    compute_architecture_family_risk(architecture_str, ...
    f_vec, params_struct, relative_sampling_iters, ...
    compute_method_flag, max_generations)

N = length(f_vec); % have a seperate function for family risk
iters = length(params_struct.z_std); % determine iterations by number of variance parameters
relative_risk = zeros(max_generations, iters);

switch compute_method_flag % how to compute
    case 'enumerate' % here compute all 2^2N genotypes for each pair of family members
        sibs_genotype = []; sibs_freqs = []; % not computed at the moment
        [x_vec p_x_vec] = initilize_x_vec_constants(N, 0, vec2row(f_vec), ...
            'enumerate', relative_sampling_iters); % take a random set of paternal genotypes
        z = genetic_architecture(x_vec, architecture_str, params_struct, 1); % get disease value
        mu = sum(p_x_vec .* z); % compute trait's mean (prevalence)
        family_tree = generate_family_tree(max_generations); [K U] = kinship_coefficient(family_tree);
        z2 = my_tensor_prod(z, z);  % z2 = repmat(z, 2^N, 1) .* reshape(repmat(z, 1, 2^N)', 2^(2*N), 1); % another tensor product for z
        
        for i=1:length(family_tree) % loop on all members
            T_marginal = joint_familial_marginal_probs(f_vec, ...
                K(i,end), reshape(U(i,end,:), 3, 1), 'genotype'); % get joint matrix
            [g_vec g_probs] = alleles_to_genotype(x_vec, p_x_vec);
            [x_vec2 p_x_vec2] = probs_tensor_product(g_vec, g_probs, g_vec, g_probs, T_marginal); % get joint 2^(2N) probabilities
            [x_vec2 p_x_vec2] = genotype_to_alleles(x_vec2, p_x_vec2);
            [x_vec2 sort_perm] = sortrows(x_vec2);
            family_risk(i) = sum(p_x_vec2(sort_perm) .* z2);
        end
        family_risk = family_risk ./ mu; % normalize by mu
        
        Risch_flag = 0;
        if( (N == 2) && Risch_flag )% one locus. Apply Risch's formula (analytic)
            p_x_times_z_vec = p_x_vec .* z;
            [mu V v_env v_gen mz_twin_risk] = ...
                compute_architecture_statistics(architecture_str, ...
                f_vec', params_struct, p_x_vec, p_x_times_z_vec, 1, 'enumerate');
            [v_marginal GRR_marginal p_z_x_marginal] = ...
                compute_architecture_statistics_marginal(architecture_str, ...
                f_vec, params_struct, p_x_vec, p_x_times_z_vec, 1, 'enumerate', mu);
            v_add = sum(V - v_marginal);
            v_gen = V - v_env;
            v_dom = v_gen - v_add;
            lambda_R_Risch = 1 + (v_add .* K(end,:) + v_dom .* U(end,:,3)) ./ mu.^2 ;
            
            family_risk_Risch = lambda_R_Risch .* mu
            family_risks_should_be_the_same = family_risk - family_risk_Risch
        end
        
    case {'analytic', 'sampling'} % for now we know how to compute relative risk only using sampling
        
        compute_old_way_sequential = 1; % just do sons, grand-sons etc. Faster but more limited
        %%        if(compute_old_way_sequential)
        x_vec = initilize_x_vec_constants(N, 0, vec2row(f_vec), ...
            'sampling', relative_sampling_iters); % take a random set of paternal genotypes
        y_vec = initilize_x_vec_constants(N, 0, vec2row(f_vec), ...
            'sampling', relative_sampling_iters); % take a random set of maternal genotypes
        
        z = zeros(relative_sampling_iters, iters); % z_child = z;
        I_z = cell(iters,1); % sets indices of genotype vector where father is diseased
        for j=1:iters % loop on different architectures
            cur_params_struct = params_struct; cur_params_struct.z_std = params_struct.z_std(j);
            if(isfield(cur_params_struct, 'min_freq'))
                cur_params_struct.min_freq = params_struct.min_freq(j);
                cur_params_struct.max_freq = params_struct.max_freq(j);
            end
            z(:,j) = genetic_architecture(x_vec, architecture_str, ...
                cur_params_struct, 1); % paternal disease value
            z(:,j) = rand(relative_sampling_iters,1) < z(:,j); % transfer probs to 0/1's
            I_z{j} = find(z(:,j)); % take only the indices where z=1 for father
        end
        for i=1:max_generations  % compute risk for descendents
            w = rand(relative_sampling_iters,N) < 2^(1-i);
            c_vec = x_vec .* w + y_vec .* (1-w); % genome of descendent: prob. 2^(1-i) from father, otherwise from ('unrelated' mother)
            for j=1:iters % loop on different architectures
                cur_params_struct = params_struct; cur_params_struct.z_std = params_struct.z_std(j);
                if(isfield(cur_params_struct, 'min_freq'))
                    cur_params_struct.min_freq = params_struct.min_freq(j);
                    cur_params_struct.max_freq = params_struct.max_freq(j);
                end
                if(~isempty(I_z{j})) % found indices where father is diseased
                    z_child = genetic_architecture(c_vec(I_z{j},:), architecture_str, ...
                        cur_params_struct, 1);
                    relative_risk(i,j) = mean(z_child);
                else % not enough information! report -1
                    relative_risk(i,j) = -1;
                end
            end
        end
        %%        else % here do a whole family tree. Slower but more general, complete
        
        num_blocks = 10; % divide to blocks to keep memory low
        lambda_mz = zeros(iters,1); lambda_dz = zeros(iters,1); mu = zeros(iters,1); 
        for cur_block = 1:num_blocks % loop on blocks
            do_block_familial = cur_block
            [family_tree family_genotype_vec] = ...
                simulate_family_genotypes(max_generations, f_vec, relative_sampling_iters / num_blocks); % simulate genotype vectors for entire family
            sibs_vec = reshape(family_genotype_vec(:,1,end-1:end) + ...
                family_genotype_vec(:,2,end-1:end), relative_sampling_iters / num_blocks, 2);
            
            [sibs_genotype sibs_freqs] = unique_with_counts(sibs_vec, 'rows');
            sibs_freqs = sibs_freqs ./ sum(sibs_freqs);
            
            family_size = length(family_tree); family_tree(end,end) = 1; % make self loop to highlight last person
            if(cur_block == 1) % initilize risk to zero
                family_risk = zeros(family_size, iters);
            end
            for i=1:iters % compute risk for last node given any other node
                run_iter_of_iters = [i iters]
                cur_params_struct = params_struct; cur_params_struct.z_std = params_struct.z_std(i);
                if(isfield(cur_params_struct, 'min_freq'))
                    cur_params_struct.min_freq = params_struct.min_freq(i);
                    cur_params_struct.max_freq = params_struct.max_freq(i);
                end
                z = zeros(relative_sampling_iters / num_blocks, family_size);
                z_binary = zeros(relative_sampling_iters / num_blocks, family_size+1);
                for j=1:family_size % loop on ancestral nodes
                    z(:,j) = genetic_architecture(family_genotype_vec(:,:,j), architecture_str, ...
                        cur_params_struct, 1);
                    z_binary(:,j) = genetic_architecture(family_genotype_vec(:,:,j), architecture_str, ...
                        cur_params_struct, 1, 0); % get zero/one sampled values
                    if(j == family_size) % repeat last one to get a MZ twin
                        z_binary(:,j+1) = genetic_architecture(family_genotype_vec(:,:,j), architecture_str, ...
                            cur_params_struct, 1, 0); % get zero/one sampled values
                    end
                end
                r_MZ = corr(z_binary(:,end), z_binary(:,end-1)); % compute correlation between monozigotic twins
                r_DZ = corr(z_binary(:,end), z_binary(:,end-2)); % compute correlation between dizogotic twins (siblings)
                lambda_mz(i) = lambda_mz(i) + sum(z_binary(:,end).*z_binary(:,end-1)); 
                lambda_dz(i) = lambda_dz(i) + sum(z_binary(:,end).*z_binary(:,end-2));
                mu(i) = mu(i) + mean(z_binary(:,end)); % update mean 

                %                h_from_twins =
                [H_from_twins(i) h_add_from_twins(i) h_liability_from_twins(i)]  = ...
                    family_segregation_to_heritability(family_tree, z_binary', 'liability');
                [H_from_twins_ADE(i) h_add_from_twins_ADE(i) h_liability_from_twins_ADE(i)]  = ...
                    family_segregation_to_heritability(family_tree, z_binary', 'liability', 'ADE');
                
                for j=1:family_size % loop on ancestral nodes
                    family_risk(j,i) = family_risk(j,i) + ...
                        mean( z(:,j) .* z(:,end) ) ./ mean(z(:,j));
                end
            end % loop on iters (different architecures)
        end % loop on blocks (to keep memory low)
        family_risk = family_risk ./ num_blocks; % normalize by # of blocks
        mu = mu ./ num_blocks;
        %%        end % New: compute sibling risk by simulating an entire family tree (should be different then risk for son!!)

        lambda_mz = lambda_mz ./ (relative_sampling_iters*mu.^2); 
        lambda_dz = lambda_dz ./ (relative_sampling_iters*mu.^2); % normalize lambda's 
        
        for i=1:iters % compute heritability from twins 
            compute_h_from_twins_i = i
            h_liability_from_twins(i) = ... % compute average of all blocks at the end
                twin_concordance_to_heritability(lambda_mz(i), lambda_dz(i), mu(i), 'ACE');
            h_liability_from_twins_ADE(i) = ... % compute average of all blocks at the end
                twin_concordance_to_heritability(lambda_mz(i), lambda_dz(i), mu(i), 'ADE');
        end
        H_from_twins = h_liability_from_twins; 
        h_add_from_twins = h_liability_from_twins;
        H_from_twins_ADE = h_liability_from_twins_ADE; 
        h_add_from_twins_ADE = h_liability_from_twins_ADE;
        
end % switch method

