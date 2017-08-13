% Test the likelihood computations given phenotype and genotype data for two-class model
%
% Input:
% N - population size
% L - # of loci
% num_individuals - sample size
% mu - mutation rate
% trait_struct - structure describing phenotype: type - quantitative or diseaes
%                                                prevalence -  prevalence (for disease trait)
% rare_cumulative_per_gene - cumulative allele freqeuncy of rare variants per gene
% theta - effective mutation rate
% full_flag - simulate and compute likelihood on full genotypes (default) or just summary statistics
% figs_dir - where to save figures
%
function test_two_class_model_likelihood(N, L, num_individuals, mu, trait_struct, ...
    rare_cumulative_per_gene, target_size_by_class_vec, theta, poisson_flag, iters, full_flag, figs_dir)

Assign24MammalsGlobalConstants;
if(~exist('full_flag', 'var') || isempty(full_flag)) % default: use entire genotypes and phenotypes (not just summary statistics)
    full_flag = 1;
end
if(~exist('poisson_flag', 'var') || isempty(poisson_flag))
    poisson_flag=0;
end
if(poisson_flag)
    expand_format_flag = 'summary';
else
    expand_format_flag = 'individual';
end

test_likelihood_sum_to_one = 0;
plot_phenotype=0;
test_s=0;
test_alpha=1;
test_beta=0;
test_alpha_beta=0; % explore alpha-beta plain
simulate_genotypes_phenotypes=0;
compute_power=0;

s_null = -0.05; % 0000000000000001; % take very low s to be close to neutral
%    s_null = -0.001; % take rather strong selection
alpha = 0.5; % 0.00000000000001; % 0.49999999999; % 0499999; % 000999999999999999; % 0.899999999999; %8333; % 0.99999999999; % 333333333; % 0.999999999999; % 0.9999999; % fraction of null variants
beta = 2; % effect size of null variants

if(test_likelihood_sum_to_one)    % New: just test that the total likelihood summing up over all genotypes X's is one
    small_n = 3; % # of individuals
    small_L=2; % # of loci
    include_phenotype=0; % full_flag=1;
    clear total_prob;
    for x2_vec = 0: (2^small_n-1)*(small_L-1) % loop for L=2
        for x_vec = 0:2^small_n-1 % test only genotype part
            if(small_L==1)
                X = my_dec2base(x_vec, 2, small_n)';
            else
                X = [my_dec2base(x_vec, 2, small_n)' my_dec2base(x2_vec, 2, small_n)'];
            end
            total_prob(x_vec+1,x2_vec+1) = exp(compute_two_class_log_likelihood( ...
                s_null, 0.99999999999999999, [], target_size_by_class_vec, N, X, [], ...
                trait_struct, [], 0, 1)); % compute only genotype part
            
            for y_vec = 0:2^small_n-1 % now sum also on y,
                Y = my_dec2base(y_vec, 2, small_n);
                total_prob_with_phenotype(x_vec+1,x2_vec+1, y_vec+1) = ...
                    exp(compute_two_class_log_likelihood( ...
                    s_null, 0.99999999999999999, beta, target_size_by_class_vec, N, X, Y, ...
                    trait_struct,  beta, 1, 1)); % include also phenotype
                total_prob_phenotype_only(x_vec+1,x2_vec+1, y_vec+1) = ...
                    exp(compute_two_class_log_likelihood( ...
                    s_null, 0.99999999999999999, beta, target_size_by_class_vec, N, X, Y, ...
                    trait_struct,  beta, -1, 1)); % compute ONLY phenotype part
            end % loop on phenotype
        end
    end
    sum_probs_should_be_one = sum(total_prob(:))
    sum_with_phenotype_should_be_one = sum(total_prob_with_phenotype(:))
    sum_only_phenotypes_should_be = (2^small_n)^2
    but_it_is = sum(total_prob_phenotype_only(:))
end % test that likelihoods sum to one


[X, y, is_null_mat, f_mat, allele_type_vec, P_poly_empirical] = simulate_two_class_genotype_phenotype(s_null, alpha, beta, ...
    rare_cumulative_per_gene, target_size_by_class_vec, ...
    N, L, num_individuals, iters, trait_type, prevalence, full_flag, poisson_flag); % generate data (genotype and phenotype). We have here only missense!!


if(full_flag)
    first_X = X(:,:,1);
else
    first_X = X(:,1);
end

s_null_vec =  -logspace(-1,-4,200); % [-0.05 -0.0000000000000001]; % % -0.01:0.0005:0; % try different values of selection
alpha_vec = 0.01:0.01:1;  % try different values of fraction of rare alleles
beta_vec = beta + [-2:0.1:2]; % [1:0.01:1.4];

if(plot_phenotype)
    for i=1:2 % Plot phenotype as function of # alleles (look at one instance)
        if(i == 1) % take all alleles
            if(full_flag)
                num_alleles_vec = sum(X(:,:,1),2);
            else
                num_alleles_vec = ... % take all only null alleles num_null_alleles_vec num_people_vec L iters] = ...
                    expand_two_class_summary_statistics(X, num_individuals, expand_format_flag);
            end
            null_str = '';
        else % only take null alleles
            if(full_flag)
                num_alleles_vec = sum(X(:,:,1) .* repmat(is_null_mat(:,1), 1, num_individuals)', 2);
            else
                [~, num_alleles_vec] = ... % take only null alleles num_null_alleles_vec num_people_vec L iters] = ...
                    expand_two_class_summary_statistics(X, num_individuals, expand_format_flag);
            end
            
            null_str = 'null';
        end
        switch expand_format_flag
            case 'individual'
                figure; plot(num_alleles_vec(:,1), y(:,1), '*');
                xlabel(['num. ' null_str ' rare alleles']); ylabel('Phenotype');
                title(['Phenotype as function of # of ' null_str ' rare alleles']);
                %    figure; boxplot(y, num_alleles_vec)
                num_y = accumarray(num_alleles_vec(:,1)+1, ones(num_individuals, 1)); % ???
                mean_y = accumarray(num_alleles_vec(:,1)+1, y(:,1)) ./ num_y;
                std_y = sqrt(accumarray(num_alleles_vec(:,1)+1, y(:,1).^2)./num_y - mean_y.^2);
                figure; errorbar(0:length(mean_y)-1, mean_y, std_y)
                xlabel(['num. ' null_str ' rare alleles']); ylabel('Phenotype');
                title(['Phenotype as function of # of ' null_str ' rare alleles']);
        end % switch
    end % loop on allele type
end % if plot_phenotype

if(test_s || test_alpha)
    ttt = cputime;
    if(test_s)
        run_s_vec = s_null_vec; run_alpha_vec = alpha; var_str = 's'; alpha_or_s = s_null; alpha_or_s_vec = s_null_vec;
    else % test_alpha
        run_s_vec = s_null; run_alpha_vec = alpha_vec; var_str = '\alpha'; alpha_or_s = alpha; alpha_or_s_vec = alpha_vec;
    end % change s
    
    %     log_like_mat_change_s = ... % compute likelihood (here vary only s)
    %         compute_two_class_log_likelihood(s_null_vec, alpha, beta, N, ...
    %         X(:,:,1), y(:,1));
    %     time_s = cputime-ttt
    %     figure;  semilogx(s_null_vec, reshape(log_like_mat_change_s, length(s_null_vec), 1), '*'); hold on;
    %     line([s_null, s_null], [min(log_like_mat_change_s(:)), max(log_like_mat_change_s(:))], 'color', 'k', 'linewidth', 2);
    %     xlabel('s^*'); ylabel('LL(s^*)');
    %     title('log-likelihood as function of s^*, the selection coefficient of null alleles');
    %     my_saveas(gcf, fullfile(figs_dir, 'likelihood', 'two_class_likelihood_vary_s'), {'epsc', 'pdf', 'jpg'});
    %     [MAX_LL MAX_I] = max(log_like_mat_change_s);
    %     max_likelihood_s_null = s_null_vec(MAX_I)
    
    
    for i=1:iters % optimize for all simulations
        run_MLE_iter = i
        if(poisson_flag)
            if(sum(allele_type_vec(:,i)) == 0)
                L_vec(i) = 0;
            else
                L_vec(i) = find(allele_type_vec(:,i) > 0, 1, 'last');
            end
            cur_is_null_vec = allele_type_vec(1:L_vec(i),i);
            cur_is_null_vec(cur_is_null_vec == STOP) = 1;
            cur_is_null_vec(cur_is_null_vec == SYNONYMOUS) = 0;
            cur_is_null_vec(cur_is_null_vec == MISSENSE) = -1;
        else
            cur_is_null_vec = s_null_mat(:,i);
        end
        [log_like_mat_change_alpha_or_s P_poly] = ... % compute likelihood (here vary only alpha)
            compute_two_class_log_likelihood(run_s_vec, run_alpha_vec, beta, target_size_by_class_vec, N, ...
            X([1:L_vec(i)  end-L_vec(i)+1:end],i)   , y(:,i), [], cur_is_null_vec, 0, full_flag, num_individuals); % don't include phenotype !!
        if(i == 1)
            Q_poly = P_poly;
        else
            field_names = fields(P_poly);
            for j=1:length(field_names)
                eval_str = ['Q_poly.' field_names{j} ' = Q_poly.' field_names{j} ' + P_poly.' field_names{j} ';']; eval(eval_str);
                if(i == iters)
                    eval_str = ['Q_poly.' field_names{j} ' = Q_poly.' field_names{j} ' ./ ' num2str(iters) ';']; eval(eval_str);
                end
            end
        end % if i==1
        
        [MAX_LL MAX_I] = max(log_like_mat_change_alpha_or_s);
        max_likelihood_alpha_or_s(i) = alpha_or_s_vec(MAX_I);
    end % loop on iters
    
    % Compare empirical probabilities to likelihood probabilities
    [~, true_I]  = min(abs(alpha_or_s_vec - alpha_or_s))
    figure; hold on;
    plot(P_poly_empirical.empirical_prob_neutral_polymorphic_in_population, P_poly.prob_neutral_allele_polymorphic_in_population, '*');
    plot(repmat(P_poly_empirical.empirical_prob_null_polymorphic_in_population, length(run_s_vec), 1), ...
        P_poly.prob_null_allele_polymorphic_in_population, 'r*');
    plot(repmat(P_poly_empirical.empirical_prob_missense_polymorphic_in_population, length(alpha_or_s_vec), 1), ...
        P_poly.prob_missense_allele_polymorphic_in_population, 'g*');
    plot(P_poly_empirical.empirical_prob_missense_polymorphic_in_population, ...
        P_poly.prob_missense_allele_polymorphic_in_population(true_I), 'mo');
    

    max_x = max(P_poly_empirical.empirical_prob_neutral_polymorphic_in_population, P_poly.prob_neutral_allele_polymorphic_in_population);
    title('Empirical vs. likelihood probabilities of being polymorphic in POPULATION');
    xlabel('Empirical'); ylabel('likelihood'); legend({'neutral', 'null', 'missense'}, 4);
    plot(linspace(0, max_x, 100), linspace(0, max_x, 100), 'k--');
    figure; hold on;
    plot(P_poly_empirical.empirical_prob_neutral_polymorphic_in_sample, P_poly.prob_neutral_allele_polymorphic_in_sample, '*');
    plot(repmat(P_poly_empirical.empirical_prob_null_polymorphic_in_sample, length(run_s_vec), 1), ...
        P_poly.prob_null_allele_polymorphic_in_sample, 'r*');
    plot(repmat(P_poly_empirical.empirical_prob_missense_polymorphic_in_sample, length(alpha_or_s_vec), 1), ...
        P_poly.prob_missense_allele_polymorphic_in_sample, 'g*');
    plot(P_poly_empirical.empirical_prob_missense_polymorphic_in_sample, ...
        P_poly.prob_missense_allele_polymorphic_in_sample(true_I), 'mo');
    
    plot(P_poly_empirical.empirical_prob_neutral_polymorphic_in_sample, ...
        P_poly_empirical.empirical_prob_neutral_polymorphic_in_population .* ...
        P_poly_empirical.numeric_prob_neutral_polymorphic_in_sample_given_population, 'o', 'markersize', 8); 
    plot(P_poly_empirical.empirical_prob_null_polymorphic_in_sample, ...
        P_poly_empirical.empirical_prob_null_polymorphic_in_population .* ...
        P_poly_empirical.numeric_prob_null_polymorphic_in_sample_given_population, 'or', 'markersize', 8); 
    
    
    
    max_x = max(P_poly_empirical.empirical_prob_neutral_polymorphic_in_sample, P_poly.prob_neutral_allele_polymorphic_in_sample);
    title('Empirical vs. likelihood probabilities of being polymorphic in SAMPLE');
    xlabel('Empirical'); ylabel('likelihood'); legend({'neutral', 'null', 'missense'}, 4);
    plot(linspace(0, max_x, 100), linspace(0, max_x, 100), 'k--');
    
    
    [~, I_sample] = min(abs(P_poly.prob_missense_allele_polymorphic_in_sample -  ...
        P_poly_empirical.empirical_prob_missense_polymorphic_in_sample));
    min_alpha_or_s1 = alpha_or_s_vec(I_sample)
    [~, I_sample] = min(abs(P_poly.prob_null_allele_polymorphic_in_sample -  ...
        P_poly_empirical.empirical_prob_null_polymorphic_in_sample));
    min_alpha_or_s2 = alpha_or_s_vec(I_sample)
    [~, I_population] = min(abs(P_poly.prob_missense_allele_polymorphic_in_population -  ...
        P_poly_empirical.empirical_prob_missense_polymorphic_in_population));
    min_alpha_or_s3 = alpha_or_s_vec(I_population)
    [~, I_population] = min(abs(P_poly.prob_null_allele_polymorphic_in_population -  ...
        P_poly_empirical.empirical_prob_null_polymorphic_in_population));
    min_alpha_or_s4 = alpha_or_s_vec(I_population)
    
    
    
    
    
    figure; hold on; plot(alpha_or_s_vec, reshape(log_like_mat_change_alpha_or_s, length(alpha_or_s_vec), 1), '*');
    line([alpha_or_s, alpha_or_s], [min(log_like_mat_change_alpha_or_s(:)), max(log_like_mat_change_alpha_or_s(:))], 'color', 'k', 'linewidth', 4);
    line([max_likelihood_alpha_or_s(end), max_likelihood_alpha_or_s(end)], [min(log_like_mat_change_alpha_or_s(:)), max(log_like_mat_change_alpha_or_s(:))], ...
        'color', 'r', 'linewidth', 4);
    xlabel(var_str); ylabel(['LL(' var_str ')']); legend('likelihood', ['true-' var_str ], ['MLE-' var_str]);
    title(['log-likelihood as function of ' var_str ', the fraction of null alleles at birth']);
    my_saveas(gcf, fullfile(figs_dir, 'likelihood', ['two_class_likelihood_vary_' var_str]), {'epsc', 'pdf', 'jpg'});
    
    figure; hold on; boxplot(max_likelihood_alpha_or_s); plot(alpha_or_s, '*g');
    ylabel(var_str); legend(['true-' var_str]); title(['MLE-' var_str ' over ' num2str(iters) ' iterations']);
    y_lim = ylim(gca); y_lim(1) = min(y_lim(1), alpha_or_s); y_lim(2) = max(y_lim(2), alpha_or_s);
    ylim(y_lim);
    title(['Mean and quantiles of estimator for ' var_str ...
        '. BIAS = ' num2str(mean(max_likelihood_alpha_or_s)-alpha_or_s) ...
        ', RMSE = ' num2str( sqrt(mean((max_likelihood_alpha_or_s-alpha_or_s).^2)) ) ]); 
    figure; hold on;
    plot(sort(max_likelihood_alpha_or_s), (1:iters)./iters); line([alpha_or_s alpha_or_s], [0 1], 'color', 'k', 'linewidth', 2);
    xlabel(var_str); legend({'MLE', ['true-' var_str]}); title(['MLE-' var_str ' over ' num2str(iters) ' iterations']);
    
    ttt = cputime - ttt
    
    
end % test alpha

if(test_beta)% Debug: Test only beta. All alleles are almost neutral, and they're all functional
    log_like_mat_change_beta = ... % compute likelihood (here vary only alpha)
        compute_two_class_log_likelihood(s_null, alpha, beta_vec, N, ...
        first_X, y(:,1), [], [], [], full_flag);
    
    figure; hold on; plot(beta_vec, reshape(log_like_mat_change_beta, length(beta_vec), 1), '*');
    line([beta, beta], [min(log_like_mat_change_beta(:))-1, max(log_like_mat_change_beta(:))+1], 'color', 'k', 'linewidth', 2);
    xlabel('\beta'); ylabel('LL(\beta)');
    title('log-likelihood as function of \beta, the effect size of null alleles');
    my_saveas(gcf, fullfile(figs_dir, 'likelihood', 'two_class_likelihood_vary_beta'), {'epsc', 'pdf', 'jpg'});
    [MAX_LL MAX_I] = max(log_like_mat_change_beta);
    max_likelihood_beta = beta_vec(MAX_I)
    %    [alpha max_likelihood_alpha; beta max_likelihood_beta; s_null max_likelihood_s_null]
end % test beta


if(test_alpha_beta)
    log_like_mat_change_alpha_beta = ... % compute likelihood (here vary only alpha)
        compute_two_class_log_likelihood(s_null, alpha_vec, beta_vec, N, ...
        first_X, y(:,1), [], [], 1, full_flag);
    figure; hold on; plot3(alpha_vec, beta_vec, log_like_mat_change_alpha_beta);
    
    
end

if(simulate_genotypes_phenotypes)
    for i=1:5 % loop on different iterations
        if(~exist('max_likelihood_s_null', 'var'))
            max_likelihood_s_null = 0;
        end
        [X_fit{i} y_fit{i} is_null_mat_fit{i} f_mat_fit{i}] = ...
            simulate_two_class_genotype_phenotype(max_likelihood_s_null, max_likelihood_alpha, max_likelihood_beta, ...
            rare_cumulative_per_gene, N, L, num_individuals, 1); % generate data
    end
    figure; hold on; plot( sort(f_mat), (1:L) ./ L, 'k', 'linewidth', 3);
    for i=1:5
        plot(sort(f_mat_fit{i}), (1:L) ./ L, [color_vec(i) '--']);
    end
    legend('Orig. Data', 'Fitted Dist'); xlabel('Derived Allele Freq.'); ylabel('Freq.');
    time_s_is = time_s
end % simulate genotypes and phenotypes with fitted parameters


if(compute_power) %    compute_two_class_log_likelihood_power()  % compute power for the two class model
    alpha = 0.333333; % set alpha
    model_params = [s_null alpha beta L]; % set parameters
    switch trait_type
        case {'disease', 'binary'}
            model_params(end+1) = prevalence;
    end
    p_two_class_mat = [];
    n_cases_vec = 1000:1000:10000; n_controls_vec = n_cases_vec; n_vec = n_cases_vec + n_controls_vec;
    iters = 5; % num. of simulations to perform
    
    [power_vec p_vals_vec stat_vec] = ...
        compute_association_power(p_two_class_mat, n_cases_vec, n_controls_vec, alpha_vec, iters, ...
        'aggregation', 'two_class_likelihood_QTL',  'population', [], model_params);
end % compute power
