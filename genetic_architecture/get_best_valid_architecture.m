% Get the best architecture from all valid ones using some criteria like
% maximizing ratio of additive and total heretability
%
% Input:
% candidate_architectures - a struct with all possible architectures
% valid_inds - which ones do we consider as valid
% params_struct - architectures parameters
% f - MAF vector
% architecture_str - string showing the type of architectures
% interesting_flag - whether or not it satisfies all criteria
%
% Output:
% good_architecture - one 'best' architecture with all its parameters
%
function good_architecture = ...
    get_best_valid_architecture(candidate_architectures, valid_inds, ...
    params_struct, f, architecture_str, interesting_flag)

AssignGeneralConstants;
AssignStatsConstants;
precision = 3;

if(isfield(candidate_architectures(1), 'mu')) % candidate architectures    
    N = size(candidate_architectures(1).v_marginal,1); % marginal was transposed
    iters = length(candidate_architectures(1).mu);
    % [dummy min_ind] = ...
    %     min(max(candidate_architectures(1).v_additive_explained(valid_inds), epsilon^2) ./ ...
    %     (max(candidate_architectures(1).v_genetic(valid_inds), epsilon^2)) .* ...
    %     max(candidate_architectures(1).lods_ratio_marginal(valid_inds,:)-1,[],2)); % combination of best r and lods-ratio
    [dummy min_ind] = ...
        min(max(candidate_architectures(1).lods_ratio_marginal(:,valid_inds)-1,[],1)); % minimize lods-ratio (this pushes freq. down ..)
    min_ind = valid_inds(min_ind);
    
    good_architecture.N = N;
    good_architecture.MAF = f(1); % assume they're all equal
    good_architecture.arch = architecture_str;
    switch architecture_str % add # of clauses and #loci in clause
        case {'CNF', 'DNF', 'sum-of-ors', 'sum-of-ands', 'and-of-sigmoids', 'or-of-sigmoids'}
            good_architecture.arch = [good_architecture.arch ...
                ['(' num2str(params_struct.num_clauses) ',' ...
                num2str(params_struct.k_in_clause) ')']];
    end
    
    good_architecture.freq = candidate_architectures(1).mu(min_ind); % mean of trait (fraction of population having it, prevalence)
    good_architecture.penet = candidate_architectures(1).penetrance(min_ind);
    good_architecture.h_add  = candidate_architectures(1).h_add(min_ind);
    
    good_architecture.h  = candidate_architectures(1).h(min_ind); % heretability of trait (fraction of variance explained by genetics)
    good_architecture.h_liability  = candidate_architectures(1).h_liability(min_ind); % heretability of trait (liability scale)
    good_architecture.r = ...
        max(candidate_architectures(1).v_additive_explained(min_ind), epsilon^2) ./ ...
        max(candidate_architectures(1).v_genetic(min_ind), epsilon^2); % ratio of additive and total heritabilities 
    [tmp_round sort_perm] =  sort(double(uint32(round(1000* ...
        candidate_architectures(1).v_marginal_explained(:,min_ind) ...
        ./  max(candidate_architectures(1).V(min_ind), epsilon^2))))./1000, 'descend'); % lame fix
    good_architecture.h_i = ...
        [num2str_delim(vec2row(tmp_round(1:1)), ', ', precision) ', ..'];  % take all marginal fractions explained
    tmp_round =  sort(double(uint32(round(1000* ...
        candidate_architectures(1).lods_ratio_marginal(:,min_ind) )))./1000, 'descend');
    good_architecture.L_i =  ...
        [num2str_delim(vec2row(tmp_round(1:1)), ', ', precision) ', ..'];
    tmp_v_pairwise = candidate_architectures(1).v_pairwise_explained(:,:,min_ind) ./ ...
        max(candidate_architectures(1).V(min_ind), epsilon^2); % lame fix. If V is small we shouldn't arrive here anyway
    tmp_v_pairwise = tmp_v_pairwise - diag(diag(tmp_v_pairwise));
    good_architecture.h_ij  = max(tmp_v_pairwise(:));
    
    
    [III JJJ] = find(tmp_v_pairwise == good_architecture.h_ij, 1);
    if(isempty(III))
        save('tmp_bad_architecture.mat', 'candidate_architectures', 'min_ind', ...
            'valid_inds', 'params_struct', 'f', 'architecture_str', ...
            'interesting_flag', 'good_architecture');
        wtf_max_not_equal_max_h_ij = good_architecture.h_ij
        total_matrix = tmp_v_pairwise
        [III JJJ] = find(tmp_v_pairwise >= good_architecture.h_ij-epsilon, 1);
    end
    max_h_ij_inds_IJ(1) = III; max_h_ij_inds_IJ(2) = JJJ; %[max_h_ij_inds_IJ(1) max_h_ij_inds_IJ(2)] = find(tmp_v_pairwise == good_architecture.h_ij, 1);
    good_architecture.h_ij = num2str_delim( round(1000*good_architecture.h_ij) ./ ...
        1000, ', ', precision); % take MAXIMAL pairwise interaction explained
    tmp_lods_pairwise = candidate_architectures(1).lods_ratio_pairwise(:,:,min_ind);
    
    tmp_lods_pairwise = tmp_lods_pairwise + 1000*eye(N);
    good_architecture.L_ij = min(tmp_lods_pairwise(:)); % minimize over diagonal ??
    %for i=1:iters % size(tmp_lods_pairwise,3)
    tmp_lods_pairwise = tmp_lods_pairwise - ...
        diag(diag(tmp_lods_pairwise)); % don't maximize over the diagonal
    %end
    good_architecture.L_ij = [num2str(max(tmp_lods_pairwise(:)), precision) ', ' ...
        num2str(good_architecture.L_ij, precision)]; % take maximum and minimum (the farthest from one)
    good_architecture.num_fields_in_table = length(fieldnames(good_architecture));  % get fields going to table only in the first time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get all fields above this point
    good_architecture.iters = length(candidate_architectures);
    good_architecture.p_z_x_marginal_full = ...
        candidate_architectures(1).p_z_x_marginal(:,:,min_ind); % Fields from here do not appear in the table
    good_architecture.p_z_x_marginal = ...
        good_architecture.p_z_x_marginal_full(sort_perm(1),:);
    good_architecture.formula = candidate_architectures(min_ind).architecture_formula; % (min_ind);
    good_architecture.interesting_flag = interesting_flag;
    good_architecture.type = architecture_to_type(architecture_str);
    good_architecture.h_ij_full = ...
        candidate_architectures(1).v_pairwise_explained(:,:,min_ind) ./ ...
        max(candidate_architectures(1).V(min_ind), epsilon^2); % lame fix again
    good_architecture.h_i_full = candidate_architectures(1).v_marginal_explained(:,min_ind) ./ ...
        max(candidate_architectures(1).V(min_ind), epsilon^2); % lame fix again
    for i=1:N % fill diagonal with marginal inheritance
        good_architecture.h_ij_full(i,i) = ...
            good_architecture.h_i_full(i); % round for better display
    end
    good_architecture.h_ij_full = ...
        round(1000.*good_architecture.h_ij_full)./1000;
    good_architecture.L_ij_full = candidate_architectures(1).lods_ratio_pairwise(:,:,min_ind);
    %         good_architecture.L_ij_full = good_architecture.L_ij_full + ...
    %             good_architecture.L_ij_full'; % no need to symmetrize
    good_architecture.L_i_full = candidate_architectures(1).lods_ratio_marginal(:,min_ind);
    for i=1:N % fill diagonal with marginal log-ratios
        good_architecture.L_ij_full(i,i) = ...
            good_architecture.L_i_full(i);
    end
    good_architecture.L_ij_full = ...
        round(1000.*good_architecture.L_ij_full)./1000; % round for better display
    good_architecture.mu_pairwise_full = candidate_architectures(1).mu_pairwise(:,:,min_ind,:);
    good_architecture.mu_pairwise = ...
        reshape(good_architecture.mu_pairwise_full(max_h_ij_inds_IJ(1), ...
        max_h_ij_inds_IJ(2),1,:), 4, 1);
    good_architecture.p_z_x_pairwise = ...
        [vec2column(good_architecture.p_z_x_marginal) .* (1-good_architecture.mu_pairwise) ...
        vec2column(good_architecture.p_z_x_marginal) .* good_architecture.mu_pairwise];
    
    good_architecture.mu_given_k_ones = candidate_architectures(1).mu_given_k_ones(:,:,min_ind);
    good_architecture.mu_given_k_ones_std = candidate_architectures(1).mu_given_k_ones_std(:,:,min_ind);
    good_architecture.relative_risk = candidate_architectures(1).relative_risk(:,min_ind);
    good_architecture.family_tree = candidate_architectures(1).family_tree; % there is only one tree (indep. of architectures)
    good_architecture.family_risk = candidate_architectures(1).family_risk(:,min_ind);
            
            
    for i=1:2
        tmp_z_x_tab = good_architecture.p_z_x_marginal_full(max_h_ij_inds_IJ(i),:);
        alpha(i) = (tmp_z_x_tab(1) + tmp_z_x_tab(2)) * tmp_z_x_tab(4) / ...
            ( (tmp_z_x_tab(3) + tmp_z_x_tab(4)) * tmp_z_x_tab(2) );
    end
    % good_architecture.p_z_x_marginal_full(max_h_ij_ind_I,4) / ...
    %     good_architecture.p_z_x_marginal_full(max_h_ij_ind_I,2);
    good_architecture.mu_pairwise_expected = zeros(4,1) + good_architecture.mu_pairwise(1); % assign beta
    good_architecture.mu_pairwise_expected([2 4]) = ...
        good_architecture.mu_pairwise_expected([2 4]) * alpha(1);% compute expected Prob(z=1|x_i,x_j) for multiplicative no-epistasis model
    good_architecture.mu_pairwise_expected([3 4]) = ...
        good_architecture.mu_pairwise_expected([3 4]) * alpha(2);% compute expected Prob(z=1|x_i,x_j) for multiplicative no-epistasis model
    
    good_architecture.params_struct = params_struct;
    good_architecture.params_struct.z_std = good_architecture.params_struct.z_std(min_ind);
    good_architecture.f_vec = f; % get frequencies
    
    good_architecture.v_additive_explained = candidate_architectures(1).v_additive_explained(min_ind);
    good_architecture.v_genetic = candidate_architectures(1).v_genetic(min_ind); % copy various variances
    good_architecture.v_pairwise_explained = candidate_architectures(1).v_pairwise_explained(:,:,min_ind);
    good_architecture.v_marginal = candidate_architectures(1).v_marginal(:,min_ind);
    good_architecture.v_marginal_explained = candidate_architectures(1).v_marginal_explained(:,min_ind);

    good_architecture.H_from_twins = candidate_architectures(1).H_from_twins(min_ind);
    good_architecture.h_add_from_twins = candidate_architectures(1).h_add_from_twins(min_ind);
    good_architecture.h_liability_from_twins = candidate_architectures(1).h_liability_from_twins(min_ind);
    good_architecture.h_liability_from_twins_ADE = candidate_architectures(1).h_liability_from_twins_ADE(min_ind);
    
    good_architecture.lambda_s = good_architecture.family_risk(end-1) / good_architecture.freq; % New: add more variable (basically mostly renaming other variables)
    good_architecture.lambda_mz = good_architecture.family_risk(end) / good_architecture.freq; % New: add more variable (basically mostly renaming other variables)
    good_architecture.GRR_marginal = good_architecture.L_i_full;
    [good_architecture.lambda_s_marginal good_architecture.lambda_s_additive ...
        good_architecture.lambda_mz_additive] = ...
                genetic_relative_risk_to_heritability(good_architecture.f_vec, ...
                good_architecture.GRR_marginal, good_architecture.freq);
    good_architecture.lambda_s_multiplicative = prod(good_architecture.lambda_s_marginal); % compute manually        
    lambda_s_and_additive_difference = good_architecture.lambda_s - good_architecture.lambda_s_additive
    good_architecture.fraction_heritability_additive = good_architecture.r;
    good_architecture.fraction_lambda_s_additive = (good_architecture.lambda_s_additive-1) ./ ...
        (good_architecture.lambda_s-1);
    good_architecture.fraction_lambda_s_additive_log = log(good_architecture.lambda_s_additive) ./ ...
        log(good_architecture.lambda_s);
    
else % here we already have good architecture
    candidate_architectures = candidate_architectures(valid_inds); 
    num_architectures = length(candidate_architectures);
    identifier_str = cell(num_architectures,1); 
    for i=1:num_architectures
       identifier_str{i} = ['N_' num2str(candidate_architectures(i).N) '_' ...
           candidate_architectures(i).arch]; 
       lods_ratio_marginal_max(i) = ...
           max(candidate_architectures(i).L_i_full-1,[],1);
    end
    [~, I J] = unique(identifier_str);
    num_unique = length(I);
    for i=1:num_unique
        cur_arch_inds = find(J == i);
        [~, min_ind] = ...
            min(lods_ratio_marginal_max(cur_arch_inds)); % minimize lods-ratio (this pushes freq. down ..)
        min_ind = cur_arch_inds(min_ind);
        good_architecture(i) = candidate_architectures(min_ind);
    end    
end % if candidate (or good) architecture



