% function test_EyreWalker
AssignGeneralConstants;
figs_dir = '../../common_disease_model/figs/EyreWalker';
new_figs_dir = fullfile(figs_dir, 'new_eric');
tau_vec=1; % [0.1 0.5 1 2 10];
sigma_epsilon_vec=[0 1]; % 0.0000000000000000000001; % noise in coupling between fitness and trait
delta=1; % multiplicative coefficient relating fitness and trait
theta_vec= [3000 300 30 3]; % Gamma scale parameter
k=0.2; % Gamma shape parameter. Value given by Eyre-Walker
num_points=5000000; % number of points to simulate
N = 10000; % effective population size



sampling_test = 0; % compute variance explained for different allele frequencies for Eyre-Walker's distribution
numeric_test = 1; % run other distributions
test_mean_s = 0; % test analytic computation of mean x for s
plot_grr_to_var=0; % simple plots relating effect size (grr) to variance explained 

if(sampling_test) % Run Eyre-Walker model using sampling
    k_vec = 0.2; % [0.1 0.2 0.5 1]; % parameter of the gamma distribution
    for link_function = {'linear'} %  'sigmoid', 'linear'}
        for sigma_epsilon = sigma_epsilon_vec % loop on noise level
            for k = k_vec % loop on gamma shape parameter k 
                run_k = k
                for tau = tau_vec
                    i_theta=1;
                    for theta = theta_vec % loop on gamma scale parameter
                        [beta s f_vec V_vec f_var_bins{i_theta} f_var_hist{i_theta}] = ...
                            EyreWalker_plot_dist(tau, sigma_epsilon, delta, ... % parameters relating s and beta
                            theta, k, [], num_points, link_function{1}, figs_dir);
                        i_theta = i_theta+1;
                    end % loop on theta
                    % New: Plot Eyre-Walker's distribution, variance explained as function of S
                    figure;
                    for i_theta=1:length(theta_vec)
                        semilogx(f_var_bins{i_theta}, smooth(f_var_hist{i_theta}, 10), ...
                            color_vec(i_theta), 'linewidth', 2); hold on;
                    end
                    xlabel('Derived Allele Frequency'); ylabel('Var. Explained Density');
                    title('Dist. of Var. Expl. as func. of allele freq. in Eyre-Walker');
                    legend_vec = [repmat('$\bar{S} = ', length(theta_vec), 1) ...
                        num2str(theta_vec') repmat('$', length(theta_vec), 1)];
                    legend_vec = mat2cellRows(legend_vec);
                    legend(legend_vec, 'interpreter', 'latex');
                    my_saveas(gcf, fullfile(figs_dir, 'EyreWalker_VarExplained_tau_1_k_0.2'), {'epsc', 'pdf'});
                end % loop on tau
                close all;
            end % loop on k (of Gamma distribution)
        end % loop on sigma
    end % loop on link function
    
end % if sampling test


if(numeric_test) % New: test the relation between selection coefficient and allele frequency
    %
    beta_bins = -2:0.05:2;
    N = 10000; % Effective population size
    x_vec = (1:2*N-1)./(2*N); % vector of allele frequencies
    dist_params_struct.sigma = 0.000001; % epsilon's variance
    dist_params_struct.tau = 1; % linear relationship
    dist_params_struct.k = 2; % parameter of Gamma distribution
    dist_params_struct.theta = 300; % parameter of Gamma distribution
    coupling_str = 'eyre_walker'; % type of coupling between effect size and selection coefficient
    % g_hist = joint_selection_effect_size_allele_freq(s_bins, beta_bins, ... % allele frequency distribution
    %     N, coupling_str, dist_params_struct ); % effective population size
    % figure; imagesc(g_hist); colorbar;
    
    
    coupling_str = 'simple_eric';
    dist_params_struct.alpha_null = 0.3333333333; % proportion of null mutations
    dist_params_struct.null_beta = 1; % st.d.s effect of null mutations
    
    s_vec = [0.0001 0.001 0.01]; num_s = length(s_vec);
    p_null_given_k_n = cell(num_s,1);
    p_null_given_less_k_n = cell(num_s,1);
    for s_ind = 1:num_s
        dist_params_struct.null_s = s_vec(s_ind); % selection coefficient of null mutations
        s_bins = (0:s_vec(s_ind):2*s_vec(s_ind)); % 0.001:0.02); % different selection coefficients
        
        [g_hist s_hist beta_hist] = ...
            joint_selection_effect_size_allele_freq(s_bins, beta_bins, ... % allele frequency distribution
            N, coupling_str, dist_params_struct); % effective population size
        
        figure; imagesc_with_labels(g_hist, beta_bins, s_bins, 0.45, [], 10); colorbar; % joint distribution of effect size and selective coefficient
        colormap ('gray'); colormap(flipud(colormap));
        xlabel('Effect size (\beta)'); ylabel('Selection coefficient (s)');
        title('Joint distribution of effect size and selection coefficient');
        
        [s_observed_bins, s_observed_hist] = observed_selection_distribution(-s_bins, s_hist, N); % Compute h(s)
        n = 100;
        p_k_n = zeros(n-1, length(s_bins));
        p_s_given_k_n = zeros(n-1, length(s_bins));
        for k=1:n-1 % Compute conditional probability given k out of n in a sample
            run_k = k
            p_k_n(k,:) = prob_observed_allele_freq(k, n, ...
                s_observed_bins, s_observed_hist, N); % probability of getting k out of N
            p_s_given_k_n(k,:) = prob_selection_given_observed_allele_freq(s_bins, k, ...
                n, s_observed_bins, s_observed_hist, N); % this doesn't involve the allele frequencies
        end
        figure; imagesc_with_labels(p_s_given_k_n, s_bins, 1:n, [], []); colorbar; % conditional distribution of effect size for each k
        colormap ('gray'); colormap(flipud(colormap));
        xlabel('Selection coefficient (s)'); ylabel('k-out-of-n');
        title('Conditional probability of selection coefficient given observing k-of-n derived alleles');
        
        
        p_joint_s_k_n = p_s_given_k_n .* p_k_n; % joint prob. of s and k
        p_null_given_k_n{s_ind} =  p_s_given_k_n(:,2) ./ sum(p_s_given_k_n(:,1:2),2);
        p_null_given_less_k_n{s_ind} = cumsum(p_joint_s_k_n(:,2)) ./ ...
            cumsum(sum(p_joint_s_k_n(:,1:2),2));
        
        
        [p_s_given_x p_x] = prob_selection_given_allele_freq(s_bins, x_vec, s_bins, s_hist, N); % compute conditional probability given frequency x in the population
        p_joint_s_x = p_s_given_x .* p_x;
        p_null_s_given_x{s_ind} = p_s_given_x(2,:) ./ sum(p_s_given_x(1:2,:));
        p_null_s_given_less_x{s_ind} = cumsum(p_joint_s_x(2,:)) ./ ...
            cumsum(sum(p_joint_s_x(1:2,:)));
        
    end % loop on selection coefficient s
    
    
    conditional_str = 'pop'; % 'sample'
    for log_flag = 0:1
        figure; % Plot conditional probability of allele being null given allele frequency
        for s_ind=1:length(s_vec)
            switch log_flag
                case 0
                    log_str = '';
                    switch conditional_str
                        case 'pop'
                            plot(x_vec, p_null_s_given_x{s_ind}, color_vec(s_ind), 'linewidth', 2);
                        case 'sample'
                            plot( (1:n-1)./n,  p_null_given_k_n{s_ind}, color_vec(s_ind), 'linewidth', 2);
                    end
                case 1
                    log_str = '_log_scale';
                    switch conditional_str
                        case 'pop'
                            semilogx(x_vec, p_null_s_given_x{s_ind}, color_vec(s_ind), 'linewidth', 2);
                        case 'sample'
                            semilogx( (1:n-1)./n,  p_null_given_k_n{s_ind}, color_vec(s_ind), 'linewidth', 2);
                    end
            end
            hold on;
        end
        xlabel('allele freq.'); ylabel('prob. null'); legend([repmat('s=', length(s_vec), 1) num2str(s_vec')]);
        title('Prob. allele is null given allele frequency');
        my_saveas(gcf, fullfile(figs_dir, ...
            ['prob_s_null_conditional_on_allele_freq' log_str]), {'epsc', 'pdf'});
        
        figure;  % Plot conditional probability of allele being null given allele frequency <= f
        for s_ind=1:length(s_vec)
            switch log_flag
                case 0
                    switch conditional_str
                        case 'pop'
                            plot(x_vec, p_null_s_given_less_x{s_ind}, color_vec(s_ind));
                        case 'sample'
                            plot( (1:n-1)./n,  p_null_given_less_k_n{s_ind}, color_vec(s_ind));
                    end
                case 1
                    switch conditional_str
                        case 'pop'
                            semilogx(x_vec, p_null_s_given_less_x{s_ind}, color_vec(s_ind), 'linewidth', 2);
                        case 'sample'
                            semilogx( (1:n-1)./n,  p_null_given_less_k_n{s_ind}, color_vec(s_ind), 'linewidth', 2);
                    end
            end
            hold on;
        end
        xlabel('allele freq.'); ylabel('prob. null'); legend([repmat('s=', length(s_vec), 1) num2str(s_vec')]);
        title('Prob. allele is null given allele frequency <= f');
        my_saveas(gcf, fullfile(figs_dir, ...
            ['prob_s_null_conditional_on_allele_freq_smaller_than_f' log_str]), {'epsc', 'pdf'});
    end % if log_flag
    
    V = effect_size_to_mean_variance_explained(beta_bins, s_bins, beta_bins, g_hist, N)
    figure; plot(beta_bins, V); title('# loci needed as function of \beta');
    xlabel('xxx'); ylabel('# loci');
    
    all_s_bins = 0.0001:0.0001:0.05; all_S_bins = 4.*N.*all_s_bins;
    beta_bins = 0:0.01:2;
    figure; hold on;
    plotxx(all_s_bins, 2 .* (1 ./ (1-exp(all_S_bins)) + 1 ./ all_S_bins), ...
        all_S_bins, 2 .* (1 ./ (1-exp(all_S_bins)) + 1 ./ all_S_bins), ...
        {'s', 'S=4Ns'}, {'Mean het., 2 x(1-x)', ''}); % x1,y1,x2,y2,xlabels,ylabels)
    title('Mean heterozigosity for allele sampled from \phi_s');
    %    plot(all_s_bins, 2 .* (1 ./ (1-exp(all_S_bins)) + 1 ./ all_S_bins), 'linewidth', 2);
    %    xlabel('s'); ylabel('Mean het., 2 x(1-x)');
    my_saveas(gcf, fullfile(figs_dir, 'mean_heterozigosity_as_function_of_s'), {'epsc', 'pdf'});
    
    
    figure; hold on;
    one_over_heterozygosity =  1 ./ (2 .* (1 ./ (1-exp(all_S_bins)) + 1 ./ all_S_bins));
    plotxx(all_s_bins, 1 ./ (2 .* (1 ./ (1-exp(all_S_bins)) + 1 ./ all_S_bins)), ...
        all_S_bins, 2 .* (1+ 2.*all_S_bins ./ 3), ...    %%%    1 ./ (2 .* (1 ./ (1-exp(all_S_bins)) + 1 ./ all_S_bins)), ...
        {'s', 'S=4Ns'}, {'One-over Mean het., 1 / 2 x(1-x)', ''}); % x1,y1,x2,y2,xlabels,ylabels)
    title('One-over Mean heterozigosity for allele sampled from \phi_s');
    hold on; plot(0, 0, 'k');
    legend({'Taylor-expansion', 'exact'}, 4);
    %    plot(all_s_bins, 2 .* (1 ./ (1-exp(all_S_bins)) + 1 ./ all_S_bins), 'linewidth', 2);
    %    xlabel('s'); ylabel('Mean het., 2 x(1-x)');
    my_saveas(gcf, fullfile(figs_dir, 'one_over_mean_heterozigosity_as_function_of_s'), {'epsc', 'pdf'});
    
    one_over_heterozygosity_different_betas = one_over_heterozygosity' * ...
        (1./beta_bins.^2);
    figure; imagesc_with_labels(log10(one_over_heterozygosity_different_betas), beta_bins, all_s_bins, ...
        30, [], 10, 0);
    xlabel('Effect size (\beta)'); ylabel('Selection coefficient (s)');
    title('Number of loci needed to achieve 100% variance explained (log10)');
    my_saveas(gcf, fullfile(figs_dir, 'number_of_loci_needed_as_function_of_s_and_beta'), {'epsc', 'pdf'});
    
    
    i=1; absorption_time_vec = []; absorption_time_vec_weighted = [];
    absorption_vec = [0.00000001/N 1/(1000*N) 1/(10*N) 1/N 10/N];
    num_limits = length(absorption_vec);
    for i=1:num_limits
        absorption_integral_limit = absorption_vec(i);
        absorption_time_vec{i} = absorption_time_by_selection(-all_s_bins, 1, N, absorption_integral_limit);
        absorption_time_vec{i} = absorption_time_vec{i} ./ absorption_time_vec{i}(1); % normalize
        
        absorption_time_vec_weighted{i} = ...
            absorption_time_by_selection(-all_s_bins, 1, N, absorption_integral_limit, [], 'freq');
        absorption_time_vec_weighted{i} = ...
            absorption_time_vec_weighted{i} ./ absorption_time_vec_weighted{i}(1); % normalize
        
        %        i = i+1;
    end
    figure; hold on;
    for i=1:num_limits
        plot(all_s_bins, absorption_time_vec{i}, color_vec(i), 'linewidth', 3);
    end
    title('Time a new allele stays in the population as a function of selection coefficient');
    xlabel('s'); ylabel('absorption-time (normalized)');
    legend_vec = [repmat('\epsilon = ', num_limits, 1) num2str(absorption_vec')];
    legend(legend_vec); ylim([0 1.01]);
    my_saveas(gcf, fullfile(figs_dir, 'mean_absorption_time_as_function_of_s'), {'epsc', 'pdf'});
    
    figure; hold on;
    for i=1:num_limits
        plot(all_s_bins, absorption_time_vec_weighted{i}, color_vec(i), 'linewidth', 3);
    end
    title('Time a new allele stays in the population as a function of selection coefficient');
    xlabel('s'); ylabel('absorption-time (normalized, weighted)');
    legend_vec = [repmat('\epsilon = ', num_limits, 1) num2str(absorption_vec')];
    legend(legend_vec); ylim([0 1.01]);
    my_saveas(gcf, fullfile(figs_dir, 'mean_absorption_time_weighted_by_allele_freq_as_function_of_s'), ...
        {'epsc', 'pdf'});
    
    
    figure;
    for i=1:num_limits
        semilogx(all_s_bins, absorption_time_vec{i}, color_vec(i), 'linewidth', 2);  hold on;
    end
    title('Time a new allele stays in the population as a function of selection coefficient (log-scale)');
    xlabel('log (s)'); ylabel('absorption-time (normalized)'); legend(legend_vec);
    my_saveas(gcf, fullfile(figs_dir, 'mean_absorption_time_as_function_of_s_log_scale'), {'epsc', 'pdf'});
    
end % if numeric test


if(test_mean_s)     % Plot just the simple relation between selection coefficient and mean allele frequency
    N = 10000; % effective population size
    %    s_vec = -linspace(0.00000001,  0.1, 10000);
    s_vec = -logspace(log10(0.0000001),  log10(0.1), 100000);
    mu_x = mean_x_given_s(s_vec, N, 0, 0.000000001);
    mu_x_var_expl = mean_x_given_s(s_vec, N, 1);
    figure;
    semilogy(mu_x, -s_vec); hold on;
    semilogy(mu_x_var_expl, -s_vec, 'r');
    xlabel('Mean Derived-Allele-Freq.'); ylabel('Selection Coefficient');
    legend({'# Alleles', 'Var. Expl.'});
    title('Mean of allele freqeuncy density');
    
    
    % Test why the hell do we get a factor of 2
    s = 0.00000001;  mu_x_var_expl_s = mean_x_given_s(-s, N, 1)
    x_vec = (1:2*N-1)./(2*N); % vector of allele frequencies
    f_vec = exp(allele_freq_spectrum(x_vec, -s, N, 0, 'log', 1));
    f_vec = normalize_hist(x_vec, f_vec);
    figure; plot(x_vec, f_vec);
    mu_x_again = mean_hist(x_vec, f_vec)
    
    
    % Now test for allele freq. (not variance explained)
    s = 0.00000001;  mu_x_expl_s = mean_x_given_s(-s, N,0)
    x_vec = (1:2*N-1)./(2*N); % vector of allele frequencies
    f_vec = exp(allele_freq_spectrum(x_vec, -s, N, 0, 'log', 0));
    f_vec = normalize_hist(x_vec, f_vec);
    figure; plot(x_vec, f_vec)
    mu_x_again = mean_hist(x_vec, f_vec)
    
end % if test mean s





