% Find maximum likelihood parameters given a set of found loci in GWAS study.
% We have a model relating selection coefficient and effect size to allele
% frequency, and infer these parameters, including variance explained, for
% each allele, from observed allele-frequency and effect size spectrum
% 
% Input:
% x_vec - minor allele frequency
% true_beta_vec - true effect size (on continuous scale)
% noisy_beta_vec - effect size actually observed (on continuous scale)
% x_grid - a grid tiling all possible allele frequencies x in [0,0.5]
% beta_discovery_grid - minimum level of beta discovered for each x in x_grid
% N - effective population size
% mu - disease prevalence
% num_cases - # cases in sample
% num_controls - # controls in sample
% alpha - pval cutoff for declaring significance
% use_power - whether to use power censoring (default is 'ON')
% s_input - use the true s (instead of fitted s) when correcting for power
% study_type - discovery/replication/combined/other
% trait_type - quantitative or binary
% sort_flag - return values sorted (1) or in original order (0, default). 
% 
% Output:
% s - selection coefficient of each locus
% a - slope relating beta to selection coefficient
% LL - loglikelihood for each value of the slope
% V - variance explained by each locus
% V_corrected  - variance explained by each locus, corrected (for finite sampling effects?)
% V_corrected_noisy  - variance explained by each locus, corrected (noisy)
% V_corrected_strings - names of each correction type
% pow_integral - denominator of power correction integral
% power_log_vec - pointwise log(power) to detect each locus
%
function [s a LL V V_corrected V_corrected_noisy V_corrected_strings ...
    pow_integral power_log_vec] = ...
    maximize_gwas_observed_allele_freq_and_effect_size_likelihood(x_vec, true_beta_vec, noisy_beta_vec, ...
    x_grid, beta_discovery_grid, ...
    N, mu, num_cases, num_controls, alpha, use_power, s_input, study_type, trait_type, ...
    sort_flag)

if(isempty(num_controls)) % for QTLs
    num_controls = zeros(size(num_cases));
end
AssignGeneralConstants;
num_loci = size(noisy_beta_vec,1); % beta could be of 'width' 2 for binary traits 
%if(~exist('fit_mode', 'var') || isempty(fit_mode))
fit_mode = 'selection';
%end
if(~exist('use_power', 'var') || isempty(use_power))
    use_power = 1;
end
if(~exist('study_type', 'var') || isempty(study_type))
    study_type = 'replication';
end
if(~exist('sort_flag', 'var') || isempty(sort_flag)) % default: don't sort !!! 
    sort_flag = 0; 
end
switch trait_type
    case {'binary', 'Binary'} % in this case the input beta vec represents both GRR and prevalence (of size 2)
        prevalence = true_beta_vec(1,2);
        true_beta_vec = true_beta_vec(:,1);
        noisy_beta_vec = noisy_beta_vec(:,1);
        num_controls_discovery = num_controls(1);
    case {'QTL', 'Quantitative'}
        prevalence = [];
        num_controls_discovery = num_controls; % empty
end

num_cases_discovery = num_cases(1);
alpha_discovery = alpha(1);
switch study_type
    case 'discovery'
        num_cases_estimated = num_cases(1);
        num_controls_estimated = num_controls(1);
        alpha_estimated = alpha(1);
    case 'replication' % take last sample size (hopefully we've got two)
        num_cases_estimated = num_cases(end);
        switch trait_type
            case {'Binary', 'binary'}
                num_controls_estimated = num_controls(end);
            case {'QTL', 'Quantitative'}
                num_controls_estimated = num_controls;
        end
        alpha_estimated = alpha(end);
    case 'combined'
        num_cases_estimated = sum(num_cases);
        num_controls_estimated = sum(num_controls);
        alpha_estimated = prod(alpha);
end


[noisy_beta_vec sort_perm] = sort(noisy_beta_vec); % Sort values by size of beta
true_beta_vec = true_beta_vec(sort_perm); % noisy_beta_vec = noisy_beta_vec(sort_perm);
x_vec = x_vec(sort_perm); % sort by effect sizes
inv_sort_perm = inv_perm(sort_perm); % compute inverse permutation to convert back 

two_side_flag = 1;
if(two_side_flag) % collapse to MAF
    x_vec = min(x_vec, 1-x_vec);
end

max_a_vec = 1 / max(abs(noisy_beta_vec(:,1)));
max_a_vec = 0.1; % min(max_a_vec, 10)
res = 0.001; a_vec = [-max_a_vec:res:max_a_vec] % allow only negative a (means negative selection) max_a_vec]; % Try different values of a. Negative a means negative selection with increased effect size
a_vec(a_vec==0) = -0.00000000000001; % cant allow exactly zero

% s_grid = -0.001:0.000005:0; % 0.005;
s_grid = -0.0001:0.00001:0; % very small grid to speed up computations % 0.005;
%s_grid = -0.05:0.0001:0; % 0.005;
s_grid(s_grid==0) = 0.00000000000001; % cant allow exactly zero
% x_grid = 0.5 .* (1:(N)) ./ N; % 0.01:0.01:0.5; % two sided (everything folded to zero) % (1:(N-1)) ./ N; % make a grid of all possible allele frequencies

switch trait_type
    case {'QTL', 'Quantitative'}
%        beta_grid = 0:0.01:max(abs(true_beta_vec)); % effect size
        fine_beta_grid = linspace(0, 5, 10000); % maximum possible beta based on ~x=0.01
    case 'Binary'
%        beta_grid = 0:0.01:max(abs(true_beta_vec)); % effect size (GRR)
        fine_beta_grid = linspace(0, 10, 10000); % maximum possible beta based on ~x=0.01
end
var_explained_grid = logspace(-10, 0, 10000); var_explained_grid(1) = 10^(-50); % ensure it's smaller than any possible var. explaied

switch trait_type % set p_mat such that variance explained is exactly var_explained_grid
    case 'Binary'
        var_explained_grid_binary = heritability_scale_change(var_explained_grid, 'binary', 0.5);
        grr_vec = heritability_to_genetic_relative_risk(var_explained_grid, 'binary', 0.5, 0.5);
        p_mat = genetic_relative_risk_to_p_z_x_marginal(0.5, grr_vec, 0.5);
        test_stat = 'chi-square-analytic';
        test_type = 'armitage';
    case {'QTL', 'Quantitative'}
        p_mat = [repmat(0.5, length(var_explained_grid), 1) ...
            vec2column(sqrt(2.*var_explained_grid)) ...
            repmat([0 1], length(var_explained_grid), 1)]; % encode QTL parameters. Assume variance one
        test_stat = 'chi-square-QTL-analytic';
        test_type = 'single-locus';
end


if(isscalar(alpha))
    alpha = repmat(alpha, length(num_cases), 1);
end
power_log_var_explained_grid = zeros(1, length(var_explained_grid));
for i=1:length(num_cases)
    power_log_var_explained_grid = power_log_var_explained_grid + ...
        log(compute_association_power(p_mat, ...
        num_cases(i), num_controls(i), alpha(i), [], test_type, test_stat)); % vector with pre-computed power for each value of variance explained V
end
allele_freq_log_mat = zeros(length(s_grid), length(x_grid));
var_explained_log_mat = zeros(length(s_grid), length(x_grid));
for i=1:length(s_grid) % compute allele frequency distribution
    run_selection_grid_i = i
    allele_freq_log_mat(i,:) = allele_freq_spectrum(x_grid, s_grid(i), ...
        N, two_side_flag, 'log');
    var_explained_log_mat(i,:) = allele_freq_spectrum(x_grid, s_grid(i), ...
        N, two_side_flag, 'log', 1); % var_explained_log_mat(i,j) is the density of variance explained  for s(i) and x(j) (log-scale)
end
allele_freq_mat = exp(allele_freq_log_mat);
var_explained_density_mat = exp(var_explained_log_mat);

power_log_MAF_mat = zeros(num_loci, length(x_grid)); % power of each locus when changing x
power_log_beta_mat = zeros(num_loci, length(fine_beta_grid)); % power of each locus when changing beta
power_mat_binary = zeros(num_loci, length(x_grid));
observed_beta_hist = zeros(num_loci, length(fine_beta_grid));
if(use_power)
    for i=1:length(noisy_beta_vec) % loop on loci
        run_beta_grid_i = i
        switch trait_type % compute beta OBSERVED variance explained for different MAF
            case {'QTL', 'Quantitative'}
                var_explained_vec = beta_to_variance_explained(noisy_beta_vec(i), x_grid, 1, 'diploid'); % take OBSERVED var explained (why??) %            var_explained_vec = 2*noisy_beta_vec(i)^2.*x_grid.*(1-x_grid); % take OBSERVED var explained (why??)
            case {'Binary', 'binary'}
                [~, var_explained_vec] = genetic_relative_risk_to_variance_explained(...
                    x_grid, noisy_beta_vec(i), mu, 'diploid');
        end
        
        [~, ~, int_inds1, int_inds2] = ...
            intervals_intersect(var_explained_vec, var_explained_vec, ...
            var_explained_grid(1:end-1), var_explained_grid(2:end), 1, 1);
        power_log_MAF_mat(i,int_inds1) = power_log_var_explained_grid(int_inds2);
        
        switch trait_type % compute MAF  OBSERVED variance explained for different beta
            case {'QTL', 'Quantitative'}
                var_explained_vec = beta_to_variance_explained(fine_beta_grid, x_vec(i), 1, 'diploid');    %2.*fine_beta_grid.^2.*x_vec(i).*(1-x_vec(i));  % take MAF and change beta
            case {'Binary', 'binary'}
                [~, var_explained_vec] = genetic_relative_risk_to_variance_explained(...
                    x_vec(i), fine_beta_grid, mu, 'diploid');
        end
        
        [temp_var_explained_vec temp_sort_perm] = sort(var_explained_vec); % sort outside intersect function
        [~, ~, int_inds1, int_inds2] = ...
            intervals_intersect((temp_var_explained_vec), (temp_var_explained_vec), ...
            var_explained_grid(1:end-1), var_explained_grid(2:end), 1, 1);
        int_inds1 = temp_sort_perm(int_inds1); %         inv_sort_perm = inv_perm(temp_sort_perm);
        power_log_beta_mat(i, int_inds1) = power_log_var_explained_grid(int_inds2);
        power_mat_binary(i, beta_discovery_grid < noisy_beta_vec(i,1)) = 1; % determine loci with power one !
        
        observed_beta_hist(i,:) = ... % o(i,j) is prob(beta_j | beta_i)
            observed_effect_size_dist([noisy_beta_vec(i,:) prevalence], x_vec(i), ...
            test_stat, alpha_estimated, fine_beta_grid, ...
            num_cases_estimated, num_controls_estimated, study_type); % Currently use replication (?) how to deal with multiple studies here?
    end % loop on snps
else % don't use power
    power_log_MAF_mat(:) = 0;
end
power_mat = exp(power_log_MAF_mat);

%    LL_grid(i,:) = log(allele_freq_spectrum_power_corrected(x_grid, s_grid(i), N, noisy_beta_vec, ...
%        mu, num_cases, num_controls, alpha));
%end

% % % % h_normalization_log_mat = zeros( length(s_grid), length(noisy_beta_vec));
% % % % for i=1:length(s_grid)
% % % %     normalization_i = i
% % % %     for j=1:length(noisy_beta_vec) % compute multiplication g(x,s) * pow(x, beta)
% % % %         if(use_power || (j==1))
% % % %             h_normalization_log_mat(i,j) = integral_hist(x_grid, ... % log
% % % %                 allele_freq_mat(i,:) .* power_mat(j,:)); % compute integral
% % % %             %        exp(allele_freq_log_mat(i,:) + power_log_MAF_mat(j,:)))); % compute integral
% % % %         else
% % % %             h_normalization_log_mat(i,j) = h_normalization_log_mat(i,1); % all are the same
% % % %         end
% % % %     end
% % % % end
% % % % h_normalization_log_mat = log(h_normalization_log_mat);
diff_x = mean(diff(x_grid)); % get grid size
h_normalization_log_mat = log(diff_x * allele_freq_mat * power_mat'); % faster. h(s, beta) is 1 / int_{x=0}^{1/2} g(x,s) * Pow(x,beta) dx

[max_LL max_ind LL power_log_vec noisy_power_log_vec] = ...
    fit_allele_freq_likelihood( ...
    s_grid, x_grid, x_vec, true_beta_vec(:,1), noisy_beta_vec(:,1), ...
    mu, h_normalization_log_mat, N, ...
    num_cases_discovery, num_controls_discovery, alpha_discovery, test_type, test_stat, ...
    two_side_flag, fit_mode, trait_type); % temp: assume everything is with discovery sample sizes

plot_figures = 0;
switch trait_type
    case {'Binary', 'binary'}
        effect_str = 'GRR';
    case {'QTL', 'Quantitative'}
        effect_str = '\beta';
end

switch fit_mode
    case 'selection'
            a = 1; % just put some dummy
        s = s_grid(max_ind);
    case 'slope'
            a = a_vec(max_ind);
        s = a .* noisy_beta_vec;
end
if(plot_figures)
    figure;
    subplot(2,2,1); hold on;
    hist(true_beta_vec, 20); xlabel(effect_str); ylabel('freq.');
    title('Effect Size Dist.');
    subplot(2,2,4); hold on;
    hist(x_vec, 20); xlabel('MAF'); ylabel('freq.');
    title('MAF Dist.');
    subplot(2,2,2); hold on;
    plot(x_vec, true_beta_vec, '.'); xlabel('MAF'); ylabel(effect_str);
    title('Effect Size vs. MAF');
    switch fit_mode
        case 'selection'
            subplot(2,2,3); hold on; plot(s_grid, LL, '.'); ...
                title('Log-Likelihood of selection coefficient');
            xlabel('selection (s)'); ylabel('LL');
            plot(s, max_LL, 'r*'); % point to the MLE
        case 'slope'
            subplot(2,2,3); hold on; plot(a_vec, LL, '.');
            title('Log-Likelihood of selection slope');
            xlabel('slope (a)'); ylabel('LL');
            plot(a, max_LL, 'r*'); % point to the MLE
    end
end % if plot_figures


use_neutral = 0;
if(exist('s_input', 'var') && ~isempty(s_input)) % use fitted or 'true' s
    use_s = s_input;
else
    if(~use_neutral)
        use_s = s;
    else % here neutral
        use_s = 0;
    end
end

switch trait_type % determine standard deviation of observed effect size. Use estimation sample size
    case {'QTL', 'Quantitative'}
        sigma_beta = 1 ./ sqrt(2 .* num_cases_estimated .* x_vec .*(1-x_vec));
    case {'Binary', 'binary'}
        sigma_beta = genetic_relative_risk_to_confidence_interval(true_beta_vec, x_vec, mu, ...
            num_cases_estimated, num_controls_estimated); % compute st.d. for GRR estimation
end
[V V_corrected V_corrected_noisy V_corrected_mat ...
    V_corrected_strings pow_integral] = ...
    correct_variance_by_power(x_vec, true_beta_vec, noisy_beta_vec, mu, ...
    sigma_beta, use_s, s_grid, x_grid, beta_discovery_grid, fine_beta_grid, observed_beta_hist, ...
    power_log_MAF_mat, power_log_beta_mat, power_log_vec, noisy_power_log_vec, ...
    allele_freq_mat, var_explained_density_mat, h_normalization_log_mat, trait_type); % compute variance explained
if(plot_figures) % plot var explained as sensitivity to s
    figure; hold on; plot(s_grid, sum(V_corrected_mat(1,:,:),3), '.');
    title('Corrected var. explained for different s estimations');
    xlabel('s'); ylabel('V^2');
    plot(use_s, sum(V_corrected(:,1)), 'r*'); % point to the MLE
end


if(~sort_flag)
    V = V(inv_sort_perm); 
    V_corrected = V_corrected(inv_sort_perm,:); 
    V_corrected_noisy = V_corrected_noisy(inv_sort_perm,:);     
    pow_integral = pow_integral(inv_sort_perm,:); 
end








% Internal Function: Find maximum likelihood parameters given a set of found loci
%
% Input:
% s_grid - grid of possible selection coefficients
% x_grid - grid of possible MAF
% x_vec - minor allele frequency
% true_beta_vec - true effect size (on continuous scale)
% noisy_beta_vec - effect size actually observed (on continuous scale)
% h_normalization_log_mat - matrix of normalization constants: h(s, beta) is 1 / int_{x=0}^{1/2} g(x,s) * Pow(x,b) dx
% N - effective population size
% num_cases - # cases in sample
% num_controls - # controls in sample
% alpha - pval cutoff for declaring significance
% test_type - type of test (marginal/epistasis etc.)
% test_stat - test statistics 
% two_side_flag - MAF or DAF
% fit_mode - slope or constant selection
%
% Output:
% max_LL - maximum likelihood for optimal s
% max_ind - index of optimal s
% LL - loglikelihood for each value of the slope
% power_log_vec - log of power to detect each locus based on true effect size
% noisy_power_log_vec - log of power to detect each locus based on observed effect size
%
function [max_LL max_ind LL power_log_vec noisy_power_log_vec] = ...
    fit_allele_freq_likelihood( ...
    s_grid, x_grid, x_vec, true_beta_vec, noisy_beta_vec, mu, ...
    h_normalization_log_mat, N, ....
    num_cases, num_controls, alpha, test_type, test_stat, ...
    two_side_flag, fit_mode, trait_type)

AssignGeneralConstants;
num_loci = length(true_beta_vec);
switch fit_mode
    case 'selection' % here fit one selection coefficient for all points (should be easier)
        LL_mat = zeros(length(s_grid), num_loci);
        
        switch trait_type
            case {'QTL', 'Quantitative'}
                p_mat = [vec2column(x_vec) vec2column(true_beta_vec) ...
                    repmat([0 1], length(x_vec), 1)]; % assume variance is one
            case {'Binary', 'binary'}
                p_mat = genetic_relative_risk_to_p_z_x_marginal(x_vec, true_beta_vec, mu);
        end
        power_log_vec = zeros(1, num_loci);
        if(isscalar(alpha))
            alpha = repmat(alpha, length(num_cases), 1);
        end
        for i=1:length(num_cases) % loop on possible many studies (discovery+replication)
            power_log_vec = power_log_vec + log(compute_association_power(p_mat, ...
                num_cases(i), num_controls(i), alpha(i), [], test_type, test_stat)); % , sampling_type); % multiply spectrum by power
        end
        
        switch trait_type
            case {'QTL', 'Quantitative'}
                p_mat = [vec2column(x_vec) vec2column(noisy_beta_vec) ...
                    repmat([0 1], length(x_vec), 1)]; % assume variance is one
            case {'Binary', 'binary'}
                p_mat = genetic_relative_risk_to_p_z_x_marginal(x_vec, noisy_beta_vec, mu);
        end
        noisy_power_log_vec = log(compute_association_power(p_mat, ...
            num_cases, num_controls, alpha, [], test_type, test_stat)); % , sampling_type); % multiply spectrum by power
        [~, ~, x_inds, x_inds2] = ...
            intervals_intersect([0 x_grid(1:end-1)+epsilon], x_grid, x_vec, x_vec, ...
            [], 0);
        missing_x = setdiff(1:length(x_vec), x_inds2)
        for j=1:length(missing_x) % complete others
            x_inds(end+1) = find(x_grid >= x_vec(missing_x(j)), 1);
            x_inds2(end+1) = missing_x(j);
        end
        %        power_log_vec(x_inds2) = power_log_MAF_mat(sub2ind(1:end, x_inds)); % copy from grid
        for i=1:length(s_grid) % simply loop over all possible s's
            run_s_i = i
            LL_mat(i,:) = vec2row(allele_freq_spectrum(x_vec, s_grid(i), N, ...
                two_side_flag, 'log')) + power_log_vec; % power not actually needed for optimization
            LL_mat(i,:) = LL_mat(i,:) - vec2row(h_normalization_log_mat(i,:));
            LL(i) = sum(LL_mat(i,:)); % sum over all points
        end
        
        % Alternative: use matlab to optimize function
        
    case 'slope' % here assume s linearly related to beta and fit one slope
        LL_mat = zeros(length(a_vec), length(noisy_beta_vec));
        for i=1:length(a_vec) % loop on slope a
            run_slope_i = i
            s_vec = abs(noisy_beta_vec) .* a_vec(i); % For QTLs assume balancing selection (strength depeneds on absolute value)
            LL_mat(i,:) = allele_freq_spectrum_power_corrected( ...
                x_vec, s_vec, N, noisy_beta_vec, ...
                mu, num_cases, num_controls, alpha,  two_side_flag, 'log', use_power); % The log-likelihood of each point
            for j=1:length(noisy_beta_vec) % normalize likelihood
                [s_min_diff s_min_ind] = min(abs(s_vec(j) - s_grid))
                LL_mat(i,j) = LL_mat(i,j) - h_normalization_log_mat(s_min_ind,j);
            end
            LL(i) = sum(LL_mat(i,:)); % sum over all points
        end
end % switch fit_mode
[max_LL max_ind] = max(LL);
