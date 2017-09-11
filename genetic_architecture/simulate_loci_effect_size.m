% Sample effect sizes at random according to a pop. genetic theory
%
% Input:
% num_loci - num. snps to simulate
% N - effective population size
% mu - prevalence
% theta - population mutation rate
% simulate_mode - 'selection' constant selection, or 'slope' linear relation between beta and selection
% s_or_a - selection coefficient s or slope for selection a
% num_cases - # cases for discovery
% num_controls - # controls for discovery
% alpha - pval cutoff
% num_points - # points to generate for each s (why not included in num. loci?)
% use_power - whether to include power censoring in simulation (default) or not
%
% Output:
% true_beta_vec - effect size on continuous scale
% f_vec - risk allele frequency
% grr_vec - genetic relative risk
% h_liab_vec - heritability (variance explained) for each locus
% discovery_indicator_vec - a binary vector saying for each locus if it was found
% discovery_power_vec - power to discover each locus according to its effect size
% observed_beta_vec - the observed (noisy) values of beta the effect size. (All values from above are for the true effect size )
% replication_beta_vec - another set of observed values for beta the effect size in an independent replication study
% 
function [ true_beta_vec f_vec grr_vec h_liab_vec ...
    discovery_indicator_vec discovery_power_vec ...
    observed_beta_vec replication_beta_vec] = ...
    simulate_loci_effect_size(num_loci, N, mu, theta, ...
    simulate_mode, s_or_a, num_cases, num_controls, alpha, num_points, use_power)

two_side_flag = 1;  scale_mode = 'log'; % mu = [];
sigma = 0.00000000000000000005; mu = 0.2; slope = -0.01; % variance, mean & slope of effect size (beta)
sampling_std = 0.000000000001; % 10 % variation of observed beta over true beta
if(~exist('use_power', 'var') || isempty(use_power))
    use_power = 1;
end
test_type = 'single-locus';
test_stat = 'chi-square-QTL-analytic';

% sigma = 0.00001; mu = 0.03; % temp debug
OLD_METHOD = 0;
if(~OLD_METHOD) % simulate based on Hartl's formula
    rand_model = 'UNIFORM';
    switch rand_model
        case 'NORMAL'
            true_beta_vec = sigma * randn(num_loci,num_points) + mu;   %First simulate effect size (continuous)
        case 'UNIFORM'
            true_beta_vec = rand(num_loci, num_points) .*  sigma + mu; % uniform distribution
    end
    
    switch simulate_mode
        case 'selection'
            s_vec = repmat(s_or_a, num_loci,1); % Try constant STRONG selection % slope .* abs(true_beta_vec); % given beta, s is a linear function
        case 'slope'
            num_loci = num_points; num_points = 1; % switch roles
            true_beta_vec = true_beta_vec';
            s_vec = s_or_a .* true_beta_vec; % different selection coefficient for each locus
    end
    f_vec = zeros(num_loci,num_points);    %    f_vec = rand(num_loci, 1); % randomize based on function
    x_grid = (1:N-1)./N;
    
    loci_power_indicator_mat = zeros(num_loci, num_points); % indicator variables
    for i=1:num_loci %   length(s_vec) % randomize from each s
        simulate_locus_i = i
        %        y_grid = exp(allele_freq_spectrum_power_corrected(x_grid, s_vec(i), ...
        %            N, true_beta_vec(i), mu, num_cases, num_controls, alpha, ...
        %            two_side_flag, 'log'));
        
        if(~use_power)
            y_grid = exp(allele_freq_spectrum(x_grid, s_vec(i), ...
                N, two_side_flag, 'log')); % ignore power (faster)
        else
            y_grid = exp(allele_freq_spectrum(x_grid, s_vec(i), ...
                N, two_side_flag, 'log')); % ignore power (faster). We deal with it later ...
            %             y_grid = exp(allele_freq_spectrum_power_corrected(x_grid, s_vec(i), ...
            %                 N, true_beta_vec(i), [], num_cases, num_controls, alpha, ...
            %                 two_side_flag, 'log', 1));
        end
        f_vec(i,:) = distrnd(x_grid, y_grid, num_points);
        %         if(use_power) % compute power and randomize detection
        %             power_inds = 1:num_points;
        %             while(~isempty(power_inds))
        %                 trait_type = 'QTL'; cur_num_points = length(power_inds);
        %                 run_power_points = cur_num_points
        %                 switch trait_type
        %                     case 'QTL'
        %                         p_mat = [vec2column(f_vec(i,power_inds)) ...
        %                             vec2column(true_beta_vec(i,power_inds)) ...
        %                             repmat([0 1], cur_num_points, 1)]; % assume variance one
        %                         test_type = 'single-locus';
        %                         test_stat = 'chi-square-QTL-analytic';
        %                     case 'binary'
        %                         p_mat = genetic_relative_risk_to_p_z_x_marginal(x, beta, mu);
        %                         test_type = 'armitage';
        %                         test_stat = 'chi-square-analytic';
        %                 end
        %                 [pow_vec p_vals_vec stat_vec non_centrality_parameter] = ...
        %                     compute_association_power(p_mat, num_cases, num_controls, ...
        %                     alpha, [], test_type, test_stat);
        %                 discovery_vec = (rand(cur_num_points,1) < vec2column(pow_vec));
        %
        %                 power_inds =  power_inds(find(~discovery_vec)); % keep not filled
        %                  f_vec(i,power_inds) = distrnd(x_grid, y_grid,length(power_inds));
        %
        %             end
        %         end
    end
    discovery_indicator_vec = ones(num_loci,num_points);
    discovery_power_vec = zeros(num_loci,num_points);
    
    
    observed_beta_vec = zeros(size(true_beta_vec)); % true_beta_vec + sampling_std .* ...
        % randn(size(true_beta_vec)) ./ sqrt(num_cases); % simulate noise on beta
    replication_beta_vec = observed_beta_vec; 
        
    % randomize observed beta values
    h_liab_vec = zeros(num_loci, num_points);
    grr_vec = zeros(num_loci, num_points);
    p_vec = QTL_params_to_p_mat(f_vec, true_beta_vec, ones(num_loci,num_points));
    for i=1:num_loci % length(s_vec) % compute parameters on a different scale
        observed_p_vec  = simulate_genotype_phenotype(p_vec(i,:), ...
            num_cases, num_controls, 1, 'QTL', 'allele'); % randomize beta according to 'true' randomization
        replication_beta_vec(i,:) = observed_p_vec(:,2); % simulate twice (one for replication)
        observed_p_vec  = simulate_genotype_phenotype(p_vec(i,:), ...
            num_cases, num_controls, 1, 'QTL', 'allele'); % randomize beta according to 'true' randomization
        observed_beta_vec(i,:) = observed_p_vec(:,2);
        
        h_liab_vec(i,:) = beta_to_variance_explained( ...
            true_beta_vec(i,:), f_vec(i,:), 1, 'diploid');
        grr_vec(i,:) = beta_to_genetic_relative_risk( ...
            true_beta_vec(i,:), f_vec(i,:), mu);
        
        
        cur_p_mat = [vec2column(f_vec(i,:)) vec2column(true_beta_vec(i,:)) ...
            repmat([0 1], num_points, 1)]; % use the TRUE beta for power here!
        cur_power  = compute_association_power(cur_p_mat, ...
            num_cases, num_controls, alpha, ...
            [], test_type, test_stat); % get power
        cur_p_mat = [vec2column(f_vec(i,:)) vec2column(observed_beta_vec(i,:)) ...
            repmat([0 1], num_points, 1)]; % use the OBSERVED beta for BINARY power here!
        cur_power_binary  = compute_association_power(cur_p_mat, ...
            num_cases, num_controls, alpha, ...
            [], test_type, test_stat, [], 1); % get binary power
        discovery_indicator_vec(i,:) = vec2row(cur_power_binary); % just copy power vec2row(rand(num_points,1) < vec2column(cur_power));
        discovery_power_vec(i,:) = vec2row(cur_power);
    end % loop on # loci
    
    %    XXX = h_liab_vec;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
else % use old method
    if(~exist('theta', 'var') || isempty(theta)) % default is neutral theory
        theta = 0;
    end
    dist_x = (1:N-1) ./ N; % possible allele-frequencies
    dist_p =  (1-dist_x).^(theta-1) ./ (dist_x); % frequency according to neutral theory
    f_vec = distrnd(dist_x, dist_p, num_loci, 1); % Sample allele frequency
    grr_vec = ones(num_loci,1)+0.01; % assume constant effect size
    %grr_vec = grr_vec + randn(num_loci,1) .* 0.002; % make gaussian fluctuation in effect size
    
    
    h_liab_vec = f_vec .* (1-f_vec); % try simple naive .. looks just like the complicated h_vec
    tol = 0.00001;
    y_vec = tol:tol:(0.25-tol);
    y_density_vec = 2 .* g_density((1 - sqrt(1-4.*y_vec)) ./ 2, N) ./ sqrt(1-4.*y_vec);
    
    num_bins = 100;
    subs_vec = ceil(f_vec .* num_bins);
    h_liab_binned_vec = accumarray(subs_vec, h_liab_vec);
    figure; bar((1:length(h_liab_binned_vec)) ./ num_bins, h_liab_binned_vec); title('Contribution to variance at each allele-frequency');
    xlabel('RAF'); ylabel('Var. Explained');
    
    figure; plot(log(sort(h_liab_vec, 'descend')));
    title('Effect sizes sorted');
    xlabel('Rank'); ylabel('Var. Explained');
end % if OLD_METHOD

function g_x = g_density(x,N)

c = 1 ./ (2.*(-log(1/N) + log(1-1/N)));
g_x = c ./ (x.*(1-x));

%function G_x = G_cumulative(x,N)
%c = 1 ./ (2*(-log(1/N) + log(1-1/N)));
%G_x = c.*(log(x) - log(1-x) + log(N) + log(1-1/N));

