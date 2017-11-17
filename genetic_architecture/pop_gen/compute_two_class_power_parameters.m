% Compute power parameters for alleles in the two class model
%
% Input:
% s_null_vec - vector of possible selection coefficients
% f_rare_vec - vector of possible cutoffs for calling a rare allele 'rare' and including it in test
% beta_vec - vector of possible effect sizes
% alpha_vec - vector of possible fractions of null (functional) alleles at birth
% full_enrichment_alpha_vec -  vector of possible fractions of null (functional) alleles at birth , allowing for only including them (100%)
% n_vec - sample size
% rare_cumulative_per_gene - ??? USED ???
% N - popualtion size (we assume equlibirum)
% prevalence - disease prevalence
% L - gene length
% mu - mutation rate
% p_val_cutoff_vec - cutoff for calling a test 'significant' in an association test
% show_s_null_ind - indices for values of s to compute power for
% c_cumulative - array of cumulative distributions of alleles frequencies as function of selection coefficient
% frac_null_by_freq_cumulative
% power_parameters_output_file - NEW! output tab-delimited table
% demographic_models_struct - NEW! WE allow here more population models !!!
% demographic_model_ind - NEW! index of the model to use !!
%
% Output:
% N_half_mat - matrix of N_50 values as function of selection s and cutoff value f for rare variants for a specific alpha at birth
% N_half_mat_only_null_enrichment - matrix of N_50 values as function of selection s and cutoff value f for rare variants for alpha=1 (i.e. we take only null alleles)
% N_half_mat_by_enrichment - matrix with N_50 as function of proprotion of nulls at birth alpha, and fraction of null kept
% power_mat - matrix with power as function of effect size beta and sample size n_vec
% new_power_mat - matrix with power as function of proprotion of nulls at birth alpha, and fraction of null kept
% n_samples_half_power_vec - vector of N_50 values as function of effect size beta
% f_rare_power_vec - vector of power as function of rare-allele frequency cutoff taken
%
%
% Todo: We need to change format: Produce a Nx i vector with:
% alpha, beta, mu, s, .., power or sample size
%
% % function [N_half_mat N_half_mat_only_null_enrichment N_half_mat_by_enrichment ...
% %     power_mat new_power_mat n_samples_half_power_vec f_rare_power_vec
function  [num_samples_struct power_struct] = compute_two_class_power_parameters( ... % new: allow only two structures: one for power and one for sample size
    equilibrium_parameters_output_file, ...
    power_parameters_output_file, demographic_model_ind)


% compute_two_class_power_parameters(two_class_stat_struct, ...
%     s_null_vec, f_rare, f_rare_vec, beta_vec, alpha_vec, full_enrichment_alpha_vec, ...
%     n_controls_vec, n_cases_vec, n_vec, rare_cumulative_per_gene, ...
%     w_x_null_mat, w_x_harmless, ...
%     N, prevalence, gene_length, mu, ...
%     p_val_cutoff, p_val_cutoff_vec, show_s_null_ind, c_cumulative, frac_null_by_freq_cumulative, ...
%     power_parameters_output_file, demographic_models_struct, demographic_model_ind) % compute paremters representing power

load(equilibrium_parameters_output_file);

power_struct = cell(3,1); num_samples_struct = cell(3,1);

if(exist(power_parameters_output_file, 'file')) % Why do we need to load this?? 
    load(power_parameters_output_file); %...
end

grr_vec =  beta_to_genetic_relative_risk(beta_vec, 0.001, prevalence);

[~, i_beta] = min(abs(grr_vec  - 5)); % 61; % choose specific beta. We choose grr to be 5
[~, i_n] = min(abs(n_vec - 50000)); %  25; % choose specific n
[~, f_ind] = min(abs(f_rare - f_rare_vec))
s_ind = show_s_null_ind(1);

demographic_models_struct.model_str{end+1} = 'equilibrium';
demographic_models_struct.num_models = length(demographic_models_struct.model_str);

%n_samples_half_power_vec = zeros(length(beta_vec), 1); % vector of sample size required for 50% power as function of beta







%N_half_mat = cell(2,1); % Compute and plot detection power
%N_half_mat_only_null_enrichment = cell(2,1);
ctr=5; % count how many power calculations we already did
for i_d = demographic_model_ind % demographic_model_ind % 1:demographic_models_struct.num_models
    switch demographic_models_struct.model_str{i_d}
        case 'equilibrium'
            rho_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:); % probability of nulls. This already depends on alpha
            num_null_alleles_cum_vec = c_cumulative{1}(show_s_null_ind,:); % cumulative num. null alleles at each frequency
            num_mixture_alleles_cum_vec = c_cumulative{2}(show_s_null_ind,:); % cumulative num. null alleles at each frequency
            all_num_null_alleles_cum_vec = c_cumulative{1}; % take all values of s (not just chosen)
            use_f_vec = f_rare_vec;
            use_s_vec = s_null_vec; use_s_null_inds = show_s_null_ind;
            total_num_alleles_per_chrom = two_class_stat_struct{1}.normalization_factor_x;
        otherwise
            demographic_models_struct.data{i_d} = load(demographic_models_struct.file_names{i_d}); % this gives both s-vec and cumulatives
            rho_vec = demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) ./ ...
                ( demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) + ...
                repmat(demographic_models_struct.data{i_d}.p_vec(end,:), 10, 1) ); % This should go DOWN!!!
            num_null_alleles_cum_vec = bsxfun(@rdivide, demographic_models_struct.data{i_d}.p_vec(end:-1:1,:),  ...
                demographic_models_struct.data{i_d}.p_vec(end:-1:1,end)); % normalize
            all_num_null_alleles_cum_vec = num_null_alleles_cum_vec;
            num_mixture_alleles_cum_vec = num_null_alleles_cum_vec .* rho_vec + ...
                bsxfun(@times,  num_null_alleles_cum_vec(1,:),  (1-rho_vec)); % Here take mixture of nulls and neutrals
            use_f_vec = demographic_models_struct.data{i_d}.x_vec;
            use_s_vec = demographic_models_struct.data{i_d}.s_vec; use_s_null_inds = 1:length(use_s_vec);
            total_num_alleles_per_chrom = demographic_models_struct.data{i_d}.num_alleles_per_chrom(end:-1:1); % get total # alleles per chrom per nuc.
    end
    
    
    run_heavy=1;
    if(run_heavy) % Very Heavy loop below !!!
        j_test_stat = ctr-1;
        for stat_used = {'chi-square-analytic', 'eric_crude_enrichment_analytic'}
            % New: vary s and beta. Compute PURE sample size when there is no dillution (we observe all null alleles, and only null alleles)
            num_samples_struct{j_test_stat}.num_samples_mat = zeros(length(s_null_vec), length(beta_vec)); % N_{1/2} as function of s and beta
            num_samples_struct{j_test_stat}.varied_parameters = {'s_null_vec', 'beta_vec'};
            num_samples_struct{j_test_stat}.alpha = 1; % we assume ideal case where we know all null alleles alpha_vec;
            num_samples_struct{j_test_stat}.p_val_cutoff = p_val_cutoff;
            num_samples_struct{j_test_stat}.f_rare = 1; % assume ideal case - no need for frequency cutoff f_rare_vec(f_ind);
            num_samples_struct{j_test_stat}.power = 0.9; % required power: 90%
            
            
            for i_s = 1:length(use_s_vec)
                %            prob_site_polymorphic_at_equilibrium = (2*N*mu) * 2 * absorption_time_by_selection(abs(s_null_vec(i_s)), 1, N, 1/(2*N), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start
                
                for j=1:length(beta_vec) % loop on beta and compute power. Which s to take here?
                    if(mod(j, 50) == 0)
                        sprintf('run_beta %ld, run_s %ld out of %ld', j, i_s, length(use_s_vec))
                    end
                    % %                     whos num_mixture_alleles_cum_vec
                    % %                     whos total_num_alleles_per_chrom
                    % %                     i_s_is = i_s
                    % %                     whos grr_vec
                    % %                     j_is = j
                    
                    % NEW: here treat grr as RR!!! (1+lambda)
                    
                    p_mat = zeros(1,4);
                    p_mat(4) = grr_vec(j) * prevalence *  all_num_null_alleles_cum_vec(i_s,end) .* total_num_alleles_per_chrom(i_s);
                    p_mat(2) = prevalence - p_mat(4);
                    p_mat(3) =  all_num_null_alleles_cum_vec(i_s,end) .* total_num_alleles_per_chrom(i_s) - p_mat(4);
                    p_mat(1) = 1-sum(p_mat);
                    
                    %                     p_mat = genetic_relative_risk_to_p_z_x_marginal( 1 .* ... % (2*N*mu) * 2 * gene_length .* ......
                    %                         all_num_null_alleles_cum_vec(i_s,end) .* total_num_alleles_per_chrom(i_s), ... %%%               c_cumulative{1}(i_s,end) * two_class_stat_struct{1}.normalization_factor_x(i_s) * (2*N*mu) * 2 * gene_length, ... % take the last one corresponding to all alleles
                    %                         grr_vec(j), prevalence); % pick a particular S. Here effect size isn't dilluted !!! (only allele frequency) % 1+beta_vec(j)*w_x_null_mat{1}(i_s,end)
                    
                    if(min(p_mat) < 0) % impossible sample size 
                        num_samples_struct{j_test_stat}.num_samples_mat(i_s,j) = -1;
                    else
                        [num_samples_struct{j_test_stat}.num_samples_mat(i_s,j)] = ...
                            compute_sample_size_from_power(p_mat, num_samples_struct{4}.power, 0.5, p_val_cutoff, ...
                            [], 'single-locus', stat_used{1}, []); % compute sample size for 90% power
                    end
                end % loop on effect size
            end % loop on selection coefficient
            j_test_stat = j_test_stat+1;
        end % loop on different test statistics
    end % if run heavy
    
    
    
    
    LOF_frac = 0.2; % we divide missense/LOF
    for i_p=1:1 % Compute N_1/2 for different values of s and f^* (this is a HEAVY COMPUTATION due to bisection)
        stat_ctr=1;
        for stat_used = {'chi-square-analytic', 'eric_crude_enrichment_analytic'}
            cur_ctr = ctr+3*(stat_ctr-1);
            p_val_cutoff = p_val_cutoff_vec(i_p);
            for k=1:3
                num_samples_struct{cur_ctr+k}.num_samples_mat = ...
                    zeros(length(use_f_vec), length(use_s_vec)); % sample size required for half power as function of allele frequency cutoff and selection coefficient
                num_samples_struct{cur_ctr+k}.varied_parameters = {'f_rare_vec', 's_null_vec'};
                num_samples_struct{cur_ctr+k}.p_val_cutoff = p_val_cutoff;
                
                num_samples_struct{cur_ctr+k}.beta = beta_vec(i_beta);
                num_samples_struct{cur_ctr+k}.power = 0.5;  % required power
                num_samples_struct{cur_ctr+k}.stat_used = stat_used{1};
            end
            num_samples_struct{cur_ctr+i_p}.alpha = alpha_vec; num_samples_struct{cur_ctr+i_p}.allele_type = 'missense-mixture';
            num_samples_struct{cur_ctr+i_p+1}.alpha = 1; num_samples_struct{cur_ctr+i_p+1}.allele_type = 'missense-null';
            num_samples_struct{cur_ctr+i_p+2}.alpha = 1; num_samples_struct{cur_ctr+i_p+2}.allele_type = 'LOF'; % enrichment / no enrichment!!!
            
            
            % New! This calculation is now for a gene.
            for i=1:length(use_f_vec) % loop on different rare-allele cutoffs
                if(mod(i, 10) == 0)
                    sprintf('comput N_50 sample_size_for_allele_freq_cutoff_ind %ld for demographic model %ld', i, i_d)
                end
                for j=1:length(show_s_null_ind) % 1:length(s_null_vec) % loop on different selection coefficients
                    p_mat = genetic_relative_risk_to_p_z_x_marginal((2*N*mu) * 2 * gene_length .* ...
                        ( num_null_alleles_cum_vec(j,i) * total_num_alleles_per_chrom(j) + ...
                        ((1-alpha_vec) / alpha_vec) * num_null_alleles_cum_vec(1,i) * total_num_alleles_per_chrom(1)  ) * (1-LOF_frac), ... % take all alleles (assuming nulls are alpha)
                        1+(grr_vec(i_beta)-1) * rho_vec(j,i), prevalence); % pick a particular S
                    num_samples_struct{cur_ctr+i_p}.num_samples_mat(i,use_s_null_inds(j)) = ...
                        compute_sample_size_from_power(p_mat, repmat(0.5, size(p_mat,1),1), 0.5, p_val_cutoff, ...
                        [], 'single-locus', stat_used{1}, [], [], [], 'analytic'); % Heaviest loop !!!  faster! analytic!!!
                    
                    % NEW!!! Compute also Eric's crude values
                    
                    % Compute optimal power assuming we fully enriched for null alleles
                    p_mat_only_null_enrichment = genetic_relative_risk_to_p_z_x_marginal( (2*N*mu) * 2 * gene_length .* ...
                        num_mixture_alleles_cum_vec(j,i) .* total_num_alleles_per_chrom(j) * (1-LOF_frac), ...
                        grr_vec(i_beta), prevalence); % pick a particular S% Compute power for full enrichment
                    num_samples_struct{cur_ctr+i_p+1}.num_samples_mat(i,use_s_null_inds(j)) = ...
                        compute_sample_size_from_power(p_mat_only_null_enrichment, repmat(0.5, length(alpha_vec), 1), 0.5, p_val_cutoff, ...
                        [], 'single-locus', stat_used{1}, [], [], [], 'analytic'); % Heaviest loop !!!  faster! analytic!!!
                    
                    % % %                     num_samples_me = ...
                    % % %                         compute_sample_size_from_power(p_mat_only_null_enrichment, repmat(0.5, length(alpha_vec), 1), 0.5, p_val_cutoff, ...
                    % % %                         [], 'single-locus', 'chi-square-analytic', [], [], [], 'analytic')
                    % % %                     num_samples_eric = ...
                    % % %                         compute_sample_size_from_power(p_mat_only_null_enrichment, repmat(0.5, length(alpha_vec), 1), 0.5, p_val_cutoff, ...
                    % % %                         [], 'single-locus', 'eric_crude_enrichment_analytic', [], [], [], 'analytic')
                    % % %                     compute_association_power(p_mat_only_null_enrichment, num_samples_eric, num_samples_eric, 0.0005, ... %   p_val_cutoff_vec(1), ...
                    % % %                     [], 'single-locus', 'chi-square-analytic', [], [], [])
                    % % %                      iters, test_type, test_stat, sampling_type, const_effect_flag, model_params, varargin)
                    
                    % Compute optimal power assuming we fully enriched for LOF alleles
                    p_mat_only_null_enrichment_LOF = genetic_relative_risk_to_p_z_x_marginal( (2*N*mu) * 2 * gene_length .* ...
                        num_mixture_alleles_cum_vec(j,i) .* total_num_alleles_per_chrom(j) * LOF_frac, ...
                        grr_vec(i_beta), prevalence); % pick a particular S% Compute power for full enrichment
                    num_samples_struct{cur_ctr+i_p+2}.num_samples_mat(i,use_s_null_inds(j)) = ...
                        compute_sample_size_from_power(p_mat_only_null_enrichment_LOF, repmat(0.5, length(alpha_vec), 1), 0.5, p_val_cutoff, ...
                        [], 'single-locus', stat_used{1}, [], [], [], 'analytic'); % Heaviest loop !!!  faster! analytic!!!
                    
                    
                end % loop on selection coefficient
            end % loop on rare allele cutoff
            stat_ctr=stat_ctr+1;
        end % loop on test statistic used
    end % loop on p-val cutoff
    run_model = demographic_models_struct.model_str{i_d}
    run_model_ind = i_d
end % loop on demographic models




power_struct{1}.power_mat = zeros(length(beta_vec), length(n_cases_vec)); % Power as function of effect size beta and sample size (num cases)
power_struct{1}.varied_parameters = {'beta_vec', 'n_vec'};
power_struct{1}.s = s_null_vec(s_ind);
power_struct{1}.p_val_cutoff = p_val_cutoff;
power_struct{1}.alpha = alpha_vec;
power_struct{1}.f_rare = f_rare_vec(f_ind);


num_samples_struct{1,1}.num_samples_mat = zeros(length(beta_vec), 1);
num_samples_struct{1,1}.varied_parameters = {'beta_vec'};
num_samples_struct{1,1}.alpha = alpha_vec;
num_samples_struct{1,1}.s = s_null_vec(1);
num_samples_struct{1,1}.p_val_cutoff = p_val_cutoff;
num_samples_struct{1,1}.f_rare = f_rare_vec(f_ind);
num_samples_struct{1,1}.power = 0.5; % required power
num_samples_struct{1,1}.p_val_cutoff = p_val_cutoff;

for i=1:length(beta_vec) % loop on beta and compute power. Which s to take here?
    if(mod(i, 20))
        run_beta = i
    end
    p_mat = genetic_relative_risk_to_p_z_x_marginal(c_cumulative{1}(s_ind,f_ind), ...
        1+beta_vec(i)*w_x_null_mat{1}(s_ind,f_ind), prevalence); % pick a particular S
    power_struct{1}.power_mat(i,:) = ... % p_vals_vec stat_vec non_centrality_parameter] = ...
        compute_association_power(p_mat, n_cases_vec, n_controls_vec, p_val_cutoff, [], ...
        'single-locus', 'chi-square-analytic', []);
    
    [num_samples_struct{1,1}.num_samples_mat(i)] = ... % non_centrality_parameter] = ...
        compute_sample_size_from_power(p_mat, 0.5, 0.5, p_val_cutoff, ...
        [], 'single-locus', 'chi-square-analytic', []);
    %        test_type, test_stat, sampling_type, const_effect_flag, model_params, varargin)
end % loop on beta vec


% Next change : modify f_rare_vec:
power_struct{2}.power_mat = zeros(length(f_rare_vec), 1);
power_struct{2}.varied_parameters = {'f_rare_vec'};
power_struct{2}.s = s_null_vec(1);
power_struct{2}.p_val_cutoff = p_val_cutoff;
power_struct{2}.alpha = alpha_vec;
power_struct{2}.n_cases = n_cases_vec(i_n);
power_struct{2}.n_controls = n_controls_vec(i_n);
power_struct{2}.beta = beta_vec(i_beta);
for i=1:length(f_rare_vec)
    p_mat2 = genetic_relative_risk_to_p_z_x_marginal(c_cumulative{1}(i), ...
        1+beta_vec(i_beta)*frac_null_by_freq_cumulative{1}(i), prevalence); % Look at beta = 4
    %    f_rare_power_vec(i) = ...
    power_struct{2}.power_mat(i) = ...
        compute_association_power(p_mat2, n_cases_vec(i_n), n_controls_vec(i_n), p_val_cutoff, [], ...
        'single-locus', 'chi-square-analytic', []);
end


% Next stage: Compute power as function of enrichment
[~, new_s_ind] = min((s_null_vec-0.001).^2); % Choose again selection coefficient
new_alpha_vec = 0.01:0.005:1; % This measures how many of the rare variants are functional
new_kept_alleles_vec = 0.01:0.01:1; % This measures how many alleles we kept in enrichment

new_frac_null_by_freq_cumulative = new_alpha_vec .* w_x_null_mat{1}(new_s_ind,f_ind) ./ ...
    (new_alpha_vec .* w_x_null_mat{1}(new_s_ind,f_ind) + (1-new_alpha_vec) .* w_x_harmless{1}(f_ind));

power_struct{3}.power_mat = zeros(length(new_alpha_vec), length(new_kept_alleles_vec)); % power as function of enrichment
power_struct{3}.varied_parameters = {'alpha_vec', 'f_rare_vec'};
power_struct{3}.s = s_null_vec(1);
power_struct{3}.p_val_cutoff = p_val_cutoff;
power_struct{3}.n_cases = n_cases_vec(i_n);
power_struct{3}.n_controls = n_controls_vec(i_n);
power_struct{3}.beta = beta_vec(i_beta);

num_samples_struct{2,1}.num_samples_mat = zeros(length(new_alpha_vec), length(new_kept_alleles_vec)); % N_{1/2} as function of enrichment
num_samples_struct{2,1}.varied_parameters = {'alpha_vec', 'f_rare_vec'};
num_samples_struct{2,1}.s = s_null_vec(1);
num_samples_struct{2,1}.p_val_cutoff = p_val_cutoff;
num_samples_struct{2,1}.beta = beta_vec(i_beta);
num_samples_struct{2,1}.power = 0.5; % required power
for i=1:length(new_alpha_vec) % loop on different alpha's (fraction of alleles which are null)
    sprintf('run_alpha_to_compute_power %ld out of %ld', i, length(new_alpha_vec))
    for j=1:length(new_kept_alleles_vec) % loop on different fraction of alleles kept
        p_mat = genetic_relative_risk_to_p_z_x_marginal(c_cumulative{1}(new_s_ind,f_ind)*new_kept_alleles_vec(j), ...
            1+beta_vec(i_beta)*new_frac_null_by_freq_cumulative(i), prevalence); % Look at beta = 4
        power_struct{3}.power_mat(i,j) = ...
            compute_association_power(p_mat, n_cases_vec(i_n), n_controls_vec(i_n), p_val_cutoff, [], ...
            'single-locus', 'chi-square-analytic', []);
        num_samples_struct{2,1}.num_samples_mat(i,j) = ...
            compute_sample_size_from_power(p_mat, 0.5, 0.5, p_val_cutoff, ...
            [], 'single-locus', 'chi-square-analytic', []);
    end % loop on fraction of alleles kept
end % loop on alpha vec





my_mkdir(dir_from_file_name(power_parameters_output_file));
if(exist(power_parameters_output_file, 'file')) % File name should now contain model information !!!
    save(power_parameters_output_file, '-append', 'num_samples_struct', 'power_struct'); %...
else
    save(power_parameters_output_file, 'num_samples_struct', 'power_struct');
end
%    'N_half_mat', 'N_half_mat_only_null_enrichment', 'N_half_mat_by_enrichment', ...
%    'power_mat', 'new_power_mat', 'n_samples_half_power_vec', 'f_rare_power_vec');



