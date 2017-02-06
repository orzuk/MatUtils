% Compute power to detect association between a locus and a trait.
% This will be a function available to user for computing power 
% There are two sampling possibilities:
% 1. Draw randomly from the population.
% 2. Draw cases and controls. Controls are assumed to be non-cases (NOT general population)
%
% Input:
% p_mat - matrix of joint probabilites of genotype and phenotype (determines effect size. We need simpler effect size) 
% n_cases_vec - vector of number of samples in test (cases or total).
%                   New:
%                   if this value is in [0,1], it denote desired power. Then we need to compute the matching sample size n!!!
% n_controls_vec - vector of number of samples in test (controls)
% alpha_vec - probability of false posities (a set of p-value cutoffs denoting first-type errors)
% iters - number of simulations used to estimate p-value when using simulations
% test_type - test performed (single locus, epsitasis, additive-model, aggretation etc.)
% test_stat - statistical test used (hypergeometric, chi-square etc.)
% sampling_type - sample from the population, or do case-control balanced sampling (default, more realistic)
% const_effect_flag - this flag (default is 'OFF') says that when computing
%                     power, we simply assume that the input effect size was the observed one.
%                     Therefore, power is BINARY - either we passed the signficance threshold or not.
% model_params - additional parameters for test
%
% Output:
% power_mat - power to detect effect (probability of rejecting null hypothesis).            Size: (#samples, #cutoffs) (second is #params)
% pvals_vec - (optional) empirical p-values we got. Used to create the power plot.          Size: 1Xiters
% stat_vec - (optional) empirical statistic-values we got. Used to create the power plot.   Size: 1Xiters
% non_centrality_parameter - parameter capturing effect size (for analytic tests).          Size: (#samples, #cutoffs) (second is #params)
%
function [power_mat, p_vals_vec, stat_vec, non_centrality_parameter] = ...
    compute_association_power(p_mat, n_cases_vec, n_controls_vec, alpha_vec, ...
    iters, test_type, test_stat, sampling_type, const_effect_flag, model_params, varargin)

AssignGeneralConstants;
AssignStatsConstants;

if(~exist('sampling_type', 'var') || isempty(sampling_type))
    if(isempty(strfind(test_stat, 'QTL')))
        sampling_type = 'case_control';
    else
        sampling_type = 'population';
    end
end
trait_struct = []; 
if(~isempty(strfind(test_stat, 'QTL')))
    trait_struct.type = 'QTL';
else
    trait_struct.type = 'binary';
end
if(~exist('const_effect_flag', 'var') || isempty(const_effect_flag))
    const_effect_flag = 0;
end
if(~exist('model_params', 'var') || isempty(model_params))
    model_params = [];
end
if(length(model_params)>=4)
    N = model_params(4); % set population size
end
% forward_power_flag = 1;
% if(max(n_cases_vec) < 1)
%     forward_power_flag = 0; % Here compute N needed to get a specified power!!!
% end

if(isempty(n_controls_vec)) % assume we've got only cases, or cases=controls
    switch trait_struct.type
        case {'binary', 'disease'}
            switch sampling_type
                case {'full-case-control', 'case-control'}
                    n_samples_vec = n_cases_vec - mod(n_cases_vec, 2); % cases=controls
                    n_cases_vec = n_samples_vec ./ 2;  n_controls_vec = n_samples_vec ./ 2;
                case 'population'
                    mu = model_params(3);
                    n_samples_vec = n_cases_vec;
                    n_cases_vec = n_samples_vec .* mu;
                    if(isempty(strfind(test_stat, 'analytic')))
                        n_cases_vec = round(n_cases_vec);
                    end
                    n_controls_vec = n_samples_vec - n_cases_vec;
            end
        case {'QTL', 'quantitative', 'continuous'}
            n_samples_vec = n_cases_vec;
    end
else
    n_samples_vec = n_cases_vec + n_controls_vec;     % take sum (proportions should be adjusted)
end
if(length(p_mat) > 8) % here no need to convert. Already 4xN
    p_vec = p_mat;
else
    if((~strcmp(test_type, 'epistasis') & ~strcmp(test_type, 'pairwise')) ...
            & (size(p_mat,2) == 4))
        p_vec = p_mat;
    else
        p_vec = vec2row(mat_into_vec(vec2row(p_mat))); % one vector of 1x4 or 1x8 probabilities
    end
end
num_cutoffs = length(alpha_vec); % number of different p-value cutoffs
if(num_cutoffs == 1) % enable multiple input matrices representing multiple effect sizes
    num_cutoffs = size(p_vec,1);
    alpha_vec = repmat(alpha_vec, 1, num_cutoffs); % alpha is a row vector !!!
end
num_sample_sizes = length(n_cases_vec); power_mat = zeros(num_sample_sizes, num_cutoffs);
% if(num_cutoffs == 1)
%     num_cutoffs = num_sample_sizes;
%     alpha_vec = repmat(alpha_vec, 1, num_cutoffs); % alpha is a row vector !!!
% end

% from here we should always have: alpha_vec is of size num_cutoffs.
%                                  p_vec of size ?? num_cutoffs?

if(length(p_vec) == 8) % check if model is with epistatic or not
    if(strcmp(test_type, 'epistasis') | strcmp(test_type, 'pairwise'))
        epistasis_flag_original = test_model_for_epistasis(p_vec, test_stat, epsilon);
        if(epistasis_flag_original)
            xxx = 999
        end
    end
end

p_vec_original = p_vec; % save original for analytic calculations. p_vec may now change
full_flag = 0;
switch sampling_type % adjust effect size
    case {'case_control', 'case_controls', 'case-control', 'case-controls'} % change p_mat to reflect case-control design
        switch trait_struct.type % case controls makes sense only for binary trait
            case {'binary', 'disease'}
                p_vec = pop_prob_to_case_control_prob(p_vec, n_cases_vec, n_controls_vec);
        end
    case {'all_population', 'population'} % just sample randomly from entire population (do nothing)
        
    case {'extreme', 'extreme-tails'} % take the tail of a continuous trait
        
    case {'full', 'full-case-control', 'full-case-controls', 'full-population'} % here sample full
        full_flag = 1;
end

if(length(p_vec) == 8) % check if indeed the model is epistatic or not
    if(strmatch(test_type, 'epistasis') | strmatch(test_type, 'pairwise'))
        epistasis_flag_case_control = test_model_for_epistasis(p_vec, test_stat, epsilon);
    end
end

% p_vec_cc = pop_prob_to_case_control_prob(p_vec, n_cases_vec, n_controls_vec);
% p_vec_cc6 = vec2row(mat_into_vec(expand_allele_table(vec_into_mat(p_vec_cc, 2)', 'additive'))); % expand to a 6x1 vec
% q_vec_cc = vec2row(mat_into_vec(expand_allele_table(vec_into_mat(p_vec, 2)', 'additive'))); % expand to a 6x1 vec
% q_vec_cc6 = pop_prob_to_case_control_prob(q_vec_cc, n_cases_vec, n_controls_vec);


if(const_effect_flag) % compute binary 'power' - did we pass threshold or not
    iters = 1;
end
num_params = max(1, size(p_mat,1));
if(length(p_mat) == 2) % special case !!! 2x2 matrix 
    num_params=1;
end
p_vals_vec = zeros(num_params*iters,num_sample_sizes); stat_vec = p_vals_vec;
analytic_flag = ~isempty(strfind(test_stat, 'analytic'));

block_size = 100; % number of iterations do do each time
if(~analytic_flag) % just run empirical tests
    [block_starts block_ends block_lens num_blocks] = divide_region_to_blocks(1, iters, block_size);
    non_centrality_parameter = []; % nothing to put here
else % no meaning for blocks
    num_blocks = 1;
end


if( (num_blocks > 1) && full_flag) % only for full flag for now (most problematic with space)
    for i=1:num_blocks % divide to blocks by iterations
        sprintf('run_block %ld = [%ld-%ld]', i, block_starts(i), block_ends(i))
        [~, p_vals_vec(block_starts(i):block_ends(i),:), stat_vec(block_starts(i):block_ends(i),:)] = ...
            compute_association_power(p_mat, n_cases_vec, n_controls_vec, alpha_vec, ...
            block_lens(i), test_type, test_stat, sampling_type, const_effect_flag, model_params);
    end
else % no full flag (summary statistics)
    switch lower(test_type)
        case {'single-locus', 'marginal', 'binary'}
            contigency_table_vec = zeros(iters*size(p_vec,1),size(p_vec,2));             %        contigency_table_vec = zeros(iters,4); % initilize counts
            trait_mode = 'allele';
        case {'additive', 'dominant', 'recessive', 'armitage', 'trend'} % count number of 'b' alleles
            contigency_table_vec = zeros(iters*size(p_vec,1),6); % initilize counts
            p_vec = expand_allele_table(p_vec, 'additive', '4XN');
            trait_mode = 'genotype'; % 3x2 tables
            %        p_vec = vec2row(mat_into_vec(expand_allele_table(vec_into_mat(p_vec, 2)', 'additive'))); % expand to a 6x1 vec
        case {'pairwise', 'epistasis'}
            contigency_table_vec = zeros(iters,8); % initilize counts
            trait_mode = 'epistasis'; % 'allele';
    end
    if(full_flag)
        contigency_table_vec = [];
    end
    diff_n_samples_vec = [n_samples_vec(1) vec2row(diff(n_samples_vec))];
    diff_n_cases_vec = [n_cases_vec(1) vec2row(diff(n_cases_vec))];
    switch lower(trait_struct.type)
        case {'binary', 'disease'}
            diff_n_controls_vec = [n_controls_vec(1) vec2row(diff(n_controls_vec))];
        case {'qtl', 'quantitative'}
            diff_n_controls_vec = zeros(size(diff_n_cases_vec));
    end
    
    s = warning('off', 'MATLAB:glmfit');
    
    non_centrality_parameter = zeros(num_sample_sizes, num_cutoffs);
    for i=1:num_sample_sizes  % loop on sample size needed
        if( (~analytic_flag) && ~ismember(test_type, {'aggregation', 'two-class-likelihood', 'two_class_likelihood'}))
            if(~const_effect_flag)
                if(full_flag)
                    cur_contigency_table_vec = zeros(iters*size(p_vec,1), ...
                        diff_n_cases_vec(i)+diff_n_controls_vec(i),4); % why 4 at the end ???
                else
                    cur_contigency_table_vec = zeros(iters*size(p_vec,1),size(p_vec,2));
                end
                for j=1:size(p_vec,1)
                    if(full_flag)
                        cur_contigency_table_vec((j-1)*iters+1:j*iters,:,:) = ... %  cur_genotype_phenotype_prod cur_genotype_sqr] = ...
                            simulate_genotype_phenotype( ...
                            p_vec(j,:), diff_n_cases_vec(i), diff_n_controls_vec(i), ...
                            iters, trait_struct.type, trait_mode, full_flag);
                    else
                        cur_contigency_table_vec((j-1)*iters+1:j*iters,:) = ... %  cur_genotype_phenotype_prod cur_genotype_sqr] = ...
                            simulate_genotype_phenotype( ...
                            p_vec(j,:), diff_n_cases_vec(i), diff_n_controls_vec(i), ...
                            iters, trait_struct.type, trait_mode, full_flag);
                    end
                end % for loop on different p_vecs (SNPs)
            else % always 'sample' the same (true) effect size
                cur_contigency_table_vec = repmat(p_vec, ... % .* diff_n_samples_vec(i)
                    iters, 1);    % no rounding here - just multiply counts
                switch trait_struct.type
                    case {'binary', 'disease'} % compute counts
                        cur_contigency_table_vec = cur_contigency_table_vec .* diff_n_samples_vec(i);
                end
            end
            if(full_flag)
                contigency_table_vec = [contigency_table_vec cur_contigency_table_vec]; % just concatenate tables
            else
                contigency_table_vec = contigency_table_vec + cur_contigency_table_vec; % add contribution
            end
            
            %        switch trait_type         % Binary trait
            %            case 'binary'
            %                cur_contigency_table_vec = mnrnd(diff_n_samples_vec(i), p_vec, iters); % randomize and add only the difference
            %            case 'QTL'
            %                beta_estimated_vec = beta_estimated_vec + cur_contigency_table_vec(1); % take mean beta
            %                f_estimated_vec = f_estimated_vec + cur_contigency_table_vec(2); % take mean MAF
            %                 [cur_samples_QTL_vec cur_samples_genotype_vec ] = ...
            %                     MixtureOfGaussiansSimulateData([(1-f_vec)^2, 2*f_vec*(1-f_vec), f_vec^2], ...
            %                     [-beta 0 beta], [1 1 1], n_samples_vec(end)); % simulate genotype (allelic) and phenotype (QTL)
            %        end
        end % if not analytic
        if(~exist('contigency_table_vec', 'var'))
            contigency_table_vec = [];
        end
        [tmp_stat_vec, tmp_pvals_vec, power_mat(i,:), non_centrality_parameter(i,:)] = ...
            internal_assoc_test(model_params, contigency_table_vec, test_type, test_stat, trait_struct.type, ...
            p_vec, p_vec_original, num_params, ...
            n_samples_vec(i), n_cases_vec(i), n_controls_vec(i), ...
            alpha_vec, i, iters, full_flag, const_effect_flag, num_sample_sizes); % compuyte test statistic
        if(~isempty(tmp_stat_vec))
            stat_vec(:,i) = tmp_stat_vec; p_vals_vec(:,i) = tmp_pvals_vec;
        end
        
        if(mod(i,10) == 0)
            computed_power_samples_i = i
            test_type_is = test_type
        end
    end % loop on sample sizes i
end % if num_blocks > 1
if(~analytic_flag) % simulated data - compute empirical p-value
    if( (size(p_vals_vec, 2) == 1) && (size(p_vals_vec, 1) < num_cutoffs) )
        p_vals_vec = vec2row(p_vals_vec);
    end
    if(size(p_vals_vec, 1) ~= iters)
        p_vals_vec = p_vals_vec';
    end
    if(size(p_vals_vec, 2) == 1) % enable multiple cutoffs for the same p-value
        p_vals_vec = repmat(p_vals_vec, 1, num_cutoffs);
        stat_vec = repmat(stat_vec, 1, num_cutoffs);
    end
    for i=1:num_sample_sizes % compute empirical power (perhaps was already computed in function 'internal_assoc_test'
        for j=1:num_cutoffs % number of different cutoffs. Why would the p-values be different each time?
            power_mat(i,j) = sum(p_vals_vec(:,i) < alpha_vec(j)) / iters; % count how many times we detected
        end
    end
end % if analytic
%end % loop on num_samples



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Temp Obsolete Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% plot_flag = 0;
% if(plot_flag)
%     switch test_stat % Plot some figures
%         case 'chi-square'
%             degs_freedom = 1;
%             figure; hold on; % Plot chi-square statistic and see why it doesn't match ...
%             [h bin_locs] = hist_density(chi_stat_vec, 100, 'g', 1, [], 0);
%             x_vec = 0:0.01:10;
%             plot(x_vec, chi2pdf(x_vec, 1), 'r');
%             legend('sampled \chi^2', ['theoretical \chi^2(' num2str(degs_freedom) ')']);
%     end
%     figure; hold on; % Plot chi-square statistic and see why it doesn't match ...
%     [h bin_locs] = hist_density(p_vals_vec, 100, 'g', 1, [], 0);
%     x_vec = 0:0.01:1;
%     plot(x_vec, ones(length(x_vec),1), 'r');
%     legend('empirical p-values', 'theoretical null p-values');
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Temp Obsolete Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



%%%%%%%%%%%% Internal Function: compute statistic and p-value according to test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
% model_params - describe the parameters of the test
% contigency_table_vec - data in either collapsed or full representation
% test_type - type of test to perform (qtl/discrete/aggregate/epistasis ...)
% test_stat - precise test statistic used
% trait_type - quantitative (QTL) or discrete (binary/disease)
% p_vec - vector with distribution parameters
% p_vec_original - ???
% num_params - number of different input distirbutions (effect-size) to compute
% n_samples_vec - vector with number of individuals
% n_cases_vec - vector with number of affected individuals (for disease traits)
% n_controls_vec - vector with number of un-affected/control individuals (for disease traits)
% alpha_vec - fraction of null (causal) mutations
% i - index of number of samples (in n_samples_vec and other vec's)
% iters - number of times to simulate in order to compute power
% full_flag - whether data (genotype+phenotypes) is given in full or collapsed summary statistics representation
%
% Output:
% stat_vec - vector of test statistics. Size: (itersXnum_params,1)
% p_vals_vec - vector of p-values. Size: (itersXnum_params,1)
% power_mat - matrix of power. Size: (1, num_params)
% non_centrality_parameter - parameter determining power. Size: (1, num_params). (For multiple cutoffs they'll all be the same)
%
function [stat_vec, p_vals_vec, power_vec, non_centrality_parameter] = ...
    internal_assoc_test(model_params, contigency_table_vec, test_type, test_stat, trait_type, ...
    p_vec, p_vec_original, num_params, ... % f_vec, grr_vec,
    n_samples, n_cases, n_controls, ...
    alpha_vec, i, iters, full_flag, const_effect_flag, num_sample_sizes)

if(~exist('full_flag', 'var') || isempty(full_flag)) % default: data is full
    full_flag = 1;
end
if(~exist('const_effect_flag', 'var') || isempty(const_effect_flag))
    const_effect_flag = 0;
end

s = warning('off', 'MATLAB:glmfit');
analytic_flag = strfind(test_stat, 'analytic');

switch strrep(lower(test_stat), '_', '-') % set parameters for test
    case {'chi-square-analytic', 'chi_square_analytic', 'chi-square', 'eric-crude-enrichment-analytic'} % New: perform analytic computation via non-centrality parameter
        [f_vec, grr_vec, mu_vec] = p_z_x_marginal_to_genetic_relative_risk(p_vec_original);
        % % % % NEW: We don't need to convert to liability scale
        % % % %         [~, h_liability] = genetic_relative_risk_to_variance_explained(f_vec, grr_vec, mu_vec, 'diploid');
        % % % %         if(h_liability > 1)
        % % % %             error('Error! GRR and f values lead to heritability > 100% !!!');
        % % % %         end
        % % % %         beta = sqrt(h_liability ./ (2.*f_vec.*(1-f_vec)));
    case {'chi-square-QTL-analytic', 'chi-square-QTL'} % QTL
        [f_vec, beta, V] = p_mat_to_QTL_params(p_vec_original);
        V_explained = min(beta_to_variance_explained(beta, f_vec, V, 'diploid'), 1); % multiply effect by two
        
    case {'epistasis-analytic', 'pairwise-epistasis-analytic', 'pairwise_epistasis-analytic'} % here test for 2x2x2 table
        
        
    case 'LRT-analytic' % This does LRT for the MLT vs. LT model analytically
        mu = p_vec(4); h_x = model_params(1:2); N = p_vec(3);
        if(N==1)
            true_model_str = 'LT';
        else
            true_model_str = 'MLT';
        end
        [LRT_mu, LRT_var] = ...
            LRT_stat_moments_MLT(2, 1, mu, h_x, true_model_str); % here h_x is on the liability scale
        
    case {'two-class-likelihood', 'two_class_likelihood', 'two_class_likelihood_QTL', 'two-class-likelihood-QTL'} % New !!! include likelihood based test for rare alleles
        s_null_vec = model_params(1);
        alpha_vec = model_params(2);
        beta_vec = model_params(3); % set parameters
        L = model_params(4); % set # loci
        N = 10000; % set effective population size
        switch trait_type
            case {'disease', 'binary'}
                prevalence = model_params(5);
            otherwise
                prevalence = [];
        end
        rare_cumulative_per_gene = 1;
        [contigency_table_vec, y, is_null_mat, f_mat] = ...
            simulate_two_class_genotype_phenotype(s_null_vec, alpha_vec, beta_vec, rare_cumulative_per_gene, ...
            N, L, n_samples, iters, trait_type, prevalence, full_flag); % Special! Prepare at once the contigency table (full_flag=0)
        %        contigency_table_vec = [contigency_table_vec' y']'; % NEW! Don't concatenate !!! % concatenate genotype and phenotype (why???)
end % switch test_stat (for setting parameters for test)


switch strrep(lower(test_stat), '_', '-') % perform test
    % NEW! Allow case-only test, where we assume that we know the frequency of controls
    
    
    case {'chi-square-analytic', 'chi_square_analytic'} % New: perform analytic computation via non-centrality parameter
        switch test_type
            case {'single-locus', 'single_locus', 'marginal', 'binary'}
                ncp_type = []; % standard test
            case {'additive', 'dominant', 'recessive', 'armitage', 'trend'}
                ncp_type = 'ott';
            case {'armitage-eliana'}
                ncp_type = 'eliana';
            case 'case-only' % here we compute the case-only test 
                %%% NEED TO FILL 
                
                
        end % switch test type
        non_centrality_parameter = genetic_relative_risk_to_non_centrality_parameter( ...
            f_vec, grr_vec, n_cases, n_controls, mu_vec, ncp_type);
        if(~const_effect_flag) % compute (approximate) probability to get significant p-value
            power_vec = 1 - ncx2cdf(chi2inv(1 - alpha_vec,1), 1, ...
                min(500, vec2row(non_centrality_parameter))); % cut ncp at 500 to save time (ncp this high gives zero anyway and power is one)
        else % compute binary value deciding if threshold was passed
            power_vec = double(chi2inv(1 - alpha_vec,1) < non_centrality_parameter+1);
        end
    case 'eric-crude-enrichment-analytic'    % NEW! crude approximation to sample size based on Eric's formula
        switch test_type
            case {'single-locus', 'single_locus', 'marginal', 'binary'} % allow only binary test
                ncp_type = [];
        end % switch test type
        p_vec_case_control = pop_prob_to_case_control_prob(p_vec_original, n_cases, n_cases); % assume balanced study
        [f_vec_case_controls, grr_vec_case_controls] = p_z_x_marginal_to_genetic_relative_risk(p_vec_case_control);
        non_centrality_parameter = sqrt(f_vec_case_controls .* n_cases) .* (grr_vec_case_controls-1);
%        f_vec_controls = p_vec_original(:,3) ./  sum(p_vec_original(:,[1 3]));
%        f_vec_cases = p_vec_original(:,4) ./  sum(p_vec_original(:,[2 4]));
        rr_vec = p_vec_original(:,4) ./ ( sum(p_vec_original(:, [2 4])) .* sum(p_vec_original(:, [3 4])) );
        power_vec = normcdf( (norminv( alpha_vec) + (rr_vec-1) .* sqrt(f_vec .* n_cases)) ./ sqrt(rr_vec) );

        
    case  'chi-square-QTL-analytic' % New: perform analytic computation via non-centrality parameter for QTL
        % [f_vec beta V] = p_mat_to_QTL_params(p_vec_original);
        %             V_explained = beta_to_variance_explained(beta, f_vec, V, 'diploid'); % multiply effect by two
        
        non_centrality_parameter = variance_explained_to_non_centrality_parameter( ...
            V_explained,  n_samples); % Temp!!!!! * 2
        if(~const_effect_flag) % compute (approximate) probability to get significant p-value
            power_vec = 1 - ncx2cdf(chi2inv(1 - alpha_vec,1), 1, ...
                min(500, non_centrality_parameter)); % truncate ncp to save time
        else % compute binary value deciding if threshold was passed
            power_vec = double(chi2inv(1 - alpha_vec,1) < non_centrality_parameter+1);
        end
    case {'probit-analytic', 'epistasis-analytic', 'pairwise-epistasis-analytic', 'pairwise_epistasis-analytic'}
        %                [f_vec beta V] = p_mat_to_QTL_params(p_vec_original);
        %                [f_vec, GRR_vec, mu] = p_z_x_marginal_to_genetic_relative_risk(p_vec_original);
        %                V_explained = min(beta_to_variance_explained(beta, f_vec, V, 'diploid'), 1); % multiply effect by two
        
        %                 [p_x_y p_x_z p_x_x_z_expected_additive disease_grr_vec] = ... % compute expected table under no-epistasis
        %                     heritability_to_p_z_x_MLT(h_x, f_vec, 1, 1, mu);
        %                p_x_x_z_expected_additive =
        if(i == 1) % need to copmute this only once !!!
            mu = model_params(3); f_vec = model_params(5);
            h_x = model_params(1); % This is the heritability explained by
            h_x_MLT = model_params(2);
            grr_vec = heritability_to_genetic_relative_risk(h_x, 'liability', f_vec, mu);
            V_explained = genetic_relative_risk_to_variance_explained(f_vec, grr_vec, mu, 'binary');
            probit_data = n_samples .* p_mat; % assume no noise
            probit_data = [probit_data(2:2:end)' probit_data(1:2:end)'+probit_data(2:2:end)'];
            x_12_probit_matrix = [0 0; 0 1; 1 0; 1 1]; % x_1, x_2, (no x_1&x_2. Additive model)
            [BETA] = glmfit( x_12_probit_matrix, probit_data, 'binomial', 'link', 'probit');
            p_x_x_z_expected_additive = normcdf(BETA(1)+ ... % first fit LT model. Then adjust to case-control probs.
                sum(repmat(BETA(2:end), 1, 4)'.*x_12_probit_matrix,2));
            p_x_x = [1-f_vec f_vec]' * [1-f_vec f_vec]; % joint probability of x_1 and x_2 (under independence - no linkage)
            p_x_x_z_expected_additive = ( mat2vec(repmat(mat2vec(p_x_x)', 2, 1)) .* ...
                mat2vec([1-p_x_x_z_expected_additive p_x_x_z_expected_additive]') )';
            
            %                 disease_grr_again = heritability_to_genetic_relative_risk(h_x, 'liability', f_vec, mu)
            [p_x_y_observed, p_x_z_observed, p_x_x_z_observed disease_grr_vec_observed] = ... % compute expected table under no-epistasis
                heritability_to_p_z_x_MLT(h_x_MLT, f_vec, 1, N, mu); % assume k=1
            if(strcmp(sampling_type, 'case-control')) % transfer both expected and observed probabilites to case-control
                p_x_x_z_expected_additive = pop_prob_to_case_control_prob(p_x_x_z_expected_additive);
                p_x_x_z_observed = pop_prob_to_case_control_prob(p_x_x_z_observed);
            end
        end % compute expected and observed tables (only once)
        
        % Sanity check: test that 2x2 table and 2x2x2 table agree
        %                 p_x_z2 = p_x_x_z_expected_additive(1:4) + p_x_x_z_expected_additive(5:8);
        %                 sanity_error_expected = (p_x_z - reshape(p_x_z2', 2, 2)') ./ p_x_z
        %                 p_x_z_observed2 = p_x_x_z_observed(1:4) + p_x_x_z_observed(5:8);
        %                 sanity_error_observed = p_x_z_observed - reshape(p_x_z_observed2', 2, 2)'
        
        non_centrality_parameter = n_samples * ...
            sum(sum((p_x_x_z_expected_additive-p_x_x_z_observed).^2 ./ ...
            p_x_x_z_expected_additive)); % compute ncp
        
        %                 non_centrality_parameter = n_samples * ...
        %                     sum(sum(p_x_x_z_expected_additive - 2.*p_x_x_z_observed + ...
        %                     p_x_x_z_observed.^2./p_x_x_z_expected_additive)); % more accurate ncp
        
        %                 n_samples_vec(1) * ...
        %         sum(sum(expected_table - 2*observed_table + observed_table.^2./expected_table))
        
        if(~const_effect_flag) % compute (approximate) probability to get significant p-value
            power_vec = 1 - ncx2cdf(chi2inv(1 - alpha_vec,1), 1, ...
                min(500, non_centrality_parameter)); % truncate ncp to save time
        else % compute binary value deciding if threshold was passed
            power_vec = double(chi2inv(1 - alpha_vec,1) < non_centrality_parameter+1);
        end
    case 'lrt-analytic' % This does LRT for the MLT vs. LT model analytically
        
        switch test_type
            case 'case-only'
                [f_vec grr_vec mu_vec] = p_z_x_marginal_to_genetic_relative_risk(p_vec_original);
                er_vec = genetic_relative_risk_to_enrichment_over_pop(grr_vec, f_vec, mu_vec);
                q_vec = er_vec.*f_vec;
                non_centrality_parameter = 2 .* n_cases * ( q_vec .* log(q_vec ./ f_vec) + (1-q_vec) .* log((1-q_vec) ./ (1-f_vec)) ) -1; % take expectation minus one!!! (expectation is NCP + deg. freedom which is one here)
                if(~const_effect_flag) % compute (approximate) probability to get significant p-value
                    power_vec = 1 - ncx2cdf(chi2inv(1 - alpha_vec,1), 1, ...
                        min(500, non_centrality_parameter)); % truncate ncp to save time
                else % compute binary value deciding if threshold was passed
                    power_vec = double(chi2inv(1 - alpha_vec,1) < non_centrality_parameter+1); % why plus one?
                end
            otherwise % OLD LRT-ratio test
                %                p_vals_vec = 1 ./ (1+exp(stat_vec)); % Bayes-factor. What's the LRT distribution? (no parameters fitted here!). Use a bayesian prior with (0.5,0.5) s
                %                null_posterior = 1 / (1+exp(LRT_mu*n_samples));
                x_alpha = log((1-alpha_vec) ./ alpha_vec); % This is the LRT for significance level alpha
                power_vec = 1 - normcdf( (x_alpha-n_samples*LRT_mu)./ ...
                    sqrt(n_samples*LRT_var) );
                
        end
                
        
    case 'gof2d-analytic' % compute gof distribution in bins analytically
        mu = p_vec(4); h_x = model_params(1:2); N = p_vec(3);
        options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance
        mu_l = fminbnd(@(x) abs(binocdf(1-1, N, x)-(1-mu)), 0, 1, options); % find mu_l that keeps the PREVALENCE
        x_mu = norminv(1-mu); x_mu_l = norminv(1-mu_l);
        %         if(N==1)
        %             true_model_str = 'LT';
        %         else
        %             true_model_str = 'MLT';
        %         end
        
        % The bins are determined by #points so change with sample size
        points_in_bin = 50; % This is the average number of points per bin
        num_bins = max(2, floor(sqrt(n_samples/points_in_bin))); % number of bins in each axis
        num_bins = min(20, num_bins); % don't allow too many bins (saves some time)
        x_grid = [-100 norminv((1:num_bins-1) ./ num_bins) 100];
        num_points_in_bin = round(n_samples / num_bins^2); % compute actual number of points
        
        z_expected_table_LT = zeros(num_bins,num_bins);
        z_expected_table_MLT = zeros(num_bins,num_bins);
        for k=1:num_bins % loop on bins for X1
            for j=1:num_bins % loop on bins for X2
                z_expected_table_LT(k,j) = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
                    z_expected_given_two_x_MLT(x1, x2, 1, h_x(1), mu, x_mu), ...% Compute expected variance
                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)) / ...
                    quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2), ...
                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)); % normalize by density
                z_expected_table_MLT(k,j) = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
                    z_expected_given_two_x_MLT(x1, x2, N, h_x(2), mu_l, x_mu_l), ...% Compute expected variance
                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)) / ...
                    quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2), ...
                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)); % normalize by density
            end
        end
        dof = num_bins^2;
        non_centrality_parameter = sum(sum( ...
            num_points_in_bin.*(z_expected_table_MLT - z_expected_table_LT).^2 ./ ...
            (z_expected_table_LT.*(1 - z_expected_table_LT)) )); % here expected are probabilities (not counts)
        power_vec = 1-ncx2cdf(chi2inv(1-alpha_vec, dof), dof, non_centrality_parameter);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Empirical test below
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'two-class-likelihood', 'two_class_likelihood', 'two_class_likelihood_QTL', 'two-class-likelihood-QTL'}  % A test for rare variants based on two-class likelihood model
        rare_cumulative_per_gene = 1; N = 10000;  % TEMP!!! need fo fill this too
        
        stat_vec = zeros(1,iters); p_vals_vec = zeros(1,iters); % initilize vecs
        for j=1:iters % loop on number of iterations (here assume num_params = 1)
            max_LL = maximize_two_class_likelihood(s_null_vec, alpha_vec, beta_vec, ...  % LL(s, alpha, beta)
                rare_cumulative_per_gene, [], N, ...
                contigency_table_vec(:,j), y(:,j), trait_type, prevalence, is_null_mat(:,j), [], full_flag);
            max_LL_no_effect = maximize_two_class_likelihood(s_null_vec, alpha_vec, 0, ...  % LL(s, alpha, 0)
                rare_cumulative_per_gene, [], N, ...
                contigency_table_vec(:,j), y(:,j), trait_type, prevalence, is_null_mat(:,j), [], full_flag);
            stat_vec(j) = max_LL - max_LL_no_effect; % take log-likelihood ratio
            dof = 1;
            p_vals_vec(j) = chi2pdf(stat_vec(j), dof);
        end
        %        power_vec = 1-ncx2cdf(chi2inv(1-alpha_vec, dof), dof, non_centrality_parameter);
        
        
        %                [p_vals_vec(:,i) stat_vec(:,i)] =
        
        % %                 % This part is not executed!
        % %             case {'QQQQTL', 'QQQchi-square-QTL'} % Perform an empirical test of QTL. (moved this to assoc_test)
        % %                 % Formula is: 1. \hat{\sbeta} = (\sum_i x_i y_i) / \sum_i x_i^2
        % %                 %             2. \Xi^2 = \beta \sum_i x_i (2 y_i - \beta x_i) assuming
        % %                 [f_vec beta V] = p_mat_to_QTL_params(p_vec_original);
        % %                 V_explained = beta_to_variance_explained(beta, f_vec, V, 'diploid'); % multiply effect by two
        % %                 stat_vec = zeros(iters,1);
        % %                 for j=1:iters % Simulate a population and their value
        % %                     cur_contigency_table_vec = simulate_genotype_phenotype( ...
        % %                         p_vec, diff_n_cases_vec, diff_n_controls_vec, trait_type);
        % %                     beta_estimate = cur_contigency_table_vec(1);
        % %                     f_vec_estimate = cur_contigency_table_vec(2);
        % %                     stat_vec(j) = beta_estimate * (2 * genotype_phenotype_prod  - ...
        % %                         beta_esimate * genotype_sqr);
        % %                     [corr_vec p_vals_vec(j)] = corr(vec2column(cur_samples_QTL_vec), ...
        % %                         vec2column(cur_samples_genotype_vec)); % compute p-val of pearson correlation
        % %                     stat_vec(j) = corr_vec.*sqrt((n_samples-2)./(1-corr_vec.^2)); % t-statistic
        % %                 end
        
    otherwise % all discrete/QTL non-analytic tests (this should also be available for QTL!!!)
        if(i == num_sample_sizes)
            figure_flag = 0; % 1 - allow figure for last iteration
        else
            figure_flag = 0;
        end
        [p_vals_vec stat_vec] = assoc_test(...
            contigency_table_vec, test_type, test_stat, ...
            repmat(n_samples, num_params, 1), ...
            full_flag, figure_flag, model_params, alpha_vec); % Perform association test
end % switch test stat used (for running test)


if(analytic_flag)
    stat_vec = []; p_vals_vec = [];  % output something
end
if(~exist('non_centrality_parameter', 'var')) % when using empirical tests, we need to compute the ncp (assume num_params = 1)
    non_centrality_parameter = mean(stat_vec);
end
if(~exist('power_vec', 'var') || isempty(power_vec)) % need to compute power mat empirically
    power_vec = mean(double(p_vals_vec < alpha_vec)); % determine if p-val smaller than cutoff (assume num_params = 1)
end

