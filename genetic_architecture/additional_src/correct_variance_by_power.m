% Add power correction for variance explained : currently do it by stripes
%
% Input:
% x_vec - minor allele frequencies 
% beta_vec - effect size (true)
% noisy_beta_vec - effect size (observed)
% s - selection coefficient
% s_grid - grid of possible selection coefficients
% x_grid - grid of possible MAF
% beta_discovery_grid - minimal beta allowing discovery for each x
% observed_beta_hist - distribution of observed (?) beta
% power_log_mat - precomputed power matrix
% power_log_vec - log of power for each locus (based on true effect size)
% noisy_power_log_vec - log of power for each locus (based on observed effect size)
% allele_freq_mat - precomputed allele density matrix
% var_explained_density_mat - precomputed variance explained density matrix
% h_normalization_log_mat - precomputed normalization matrix, h(s, beta) is 1 / int_{x=0}^{1/2} g(x,s) * Pow(x,b) dx
%
% Output:
% V - variance of observed loci
% V_corrected - corrected variance with hidden loci
% V_corrected_noisy - ?? 
% V_corrected_mat - matrix of all corrections for different s values
% V_corrected_strings - names of different corrections
% pow_integral_ratio - fraction to multiply by for power correction
%
function [V V_corrected V_corrected_noisy V_corrected_mat V_corrected_strings ...
    pow_integral_ratio ] = ...
    correct_variance_by_power(x_vec, beta_vec, noisy_beta_vec, mu, ...
    sigma_beta, s, s_grid, x_grid, beta_discovery_grid, beta_grid, observed_beta_hist, ...
    power_log_mat, power_log_beta_mat, power_log_vec, noisy_power_log_vec, ...
    allele_freq_mat, var_explained_density_mat, h_normalization_log_mat, trait_type)


AssignGeneticArchitectureConstants;

%correction_mode = 'noisy'; % NEW! add variation in beta
%single_mean_beta = mean(beta_vec); % apply power correction using this beta for all loci

num_loci = length(x_vec);
a_limit_vec = [0 0.01 0 0.01 0 0 0.01]; % last one is with s equals zero !!! and integral goes only to x !!!
b_limit_vec = [vec2column(x_vec) ... % vec2column(x_vec) ...
    repmat(0.5, length(beta_vec), 4) vec2column(x_vec) vec2column(x_vec)];
num_corrections = length(a_limit_vec); % number of different corrections
V_corrected_strings = cell(num_corrections, 1);
switch trait_type
    case {'QTL', 'Quantitative'}
        V = vec2column(beta_to_variance_explained(beta_vec, x_vec, 1, 'diploid')); % true variance
        V_noisy = vec2column(beta_to_variance_explained(noisy_beta_vec, x_vec, 1, 'diploid')); % observed variance

    case {'Binary', 'binary'}
        [~, V] = genetic_relative_risk_to_variance_explained(x_vec, beta_vec, mu, 'diploid'); % true variance for GRR
        [~, V_noisy] = genetic_relative_risk_to_variance_explained(x_vec, noisy_beta_vec, mu, 'diploid'); % observed variance for GRR
        V = vec2column(V); V_noisy = vec2column(V_noisy); 
end
pow_integral_ratio = zeros(num_loci, length(s_grid)); % num_corrections); % pow_integral_ratio(i) = int_{x=0}^{1/2} g(x,s) * Pow(x,beta) * x(1-x) dx / int_{x=0}^{1/2} g(x,s) * x(1-x) dx
pow_integral_ratio_noisy = pow_integral_ratio;
pow_integral_only = zeros(num_loci,num_corrections);
pow_integral_one_minus_mat = zeros(num_loci, num_corrections, length(s_grid)); % store the ... % pow_integral_one_minus_mat(i) = int_{x=a}^{b} g(x,s) * (1-Pow(x,beta)) * x(1-x) dx / int_{x=0}^{1/2} g(x,s) * x(1-x) dx
pow_integral_one_minus_mat_noisy = pow_integral_one_minus_mat;
pow_integral_denominator = zeros(length(s_grid),1);
s_ind = zeros(num_corrections,1);
for j=1:num_corrections % loop on integral limits
    run_integral_limit = j
    if(ismember(j, [1 6 7]))
        V_corrected_strings{j} = ['Pop. Gen. [' num2str(a_limit_vec(j), 2) ', f]'];
    else
        V_corrected_strings{j} = ['Pop. Gen. [' num2str(a_limit_vec(j), 2) ...
            ', ' num2str(b_limit_vec(1,j), 2) ']'];
    end
    if(j < 4) % num_corrections)
        V_corrected_strings{j} = [V_corrected_strings{j} ', fitted s'];
        [~, s_ind(j)] = min(abs(s_grid == s));
    else % set neutral s
        V_corrected_strings{j} = [V_corrected_strings{j} ', s=0'];
        [~, s_ind(j)] = min(abs(s_grid)); % take s=0 (neutral) for last two
    end
    a_ind = find(x_grid > a_limit_vec(j), 1);
    
    for ss_ind = 1:length(s_grid) % New: loop on s to determine sensitivity to fitting its value
        pow_integral_denominator(ss_ind) = integral_hist(x_grid, ...
            var_explained_density_mat(ss_ind,:));
    end
    for i=1:length(beta_vec) % loop on effect size
        run_beta_ind = i
        b_ind = find(x_grid <= b_limit_vec(i,j), 1, 'last');
        %        convolved_power_and_beta_vec = exp(power_log_mat(i,:)) .* ...            % Perform convulution of noisy beta and other elements
        for ss_ind = 1:length(s_grid) % New: loop on s to determine sensitivity to fitting its value
            %            run_s_grid = ss_ind
            debug_plots = 0;
            if(debug_plots)
                if((i == 1) && (j == 3)) % temp debug plots
                    normalization_const = integral_hist(x_grid, allele_freq_mat(s_ind(j),:) .* x_grid .* (1-x_grid))
                    var_exp_grid = normalize_hist(x_grid, allele_freq_mat(s_ind(j),:) .* x_grid .* (1-x_grid));
                    var_exp_power_grid = allele_freq_mat(s_ind(j),:) .* x_grid .* (1-x_grid) .* exp(power_log_mat(i,:)) ./ normalization_const;
                    var_exp_power_grid_normalized = normalize_hist(x_grid, var_exp_power_grid);
                    
                    normalization_const = integral_hist(x_grid, allele_freq_mat(s_ind(j),:))
                    allele_freq_grid = normalize_hist(x_grid, allele_freq_mat(s_ind(j),:));
                    allele_freq_power_grid = allele_freq_mat(s_ind(j),:)  .* exp(power_log_mat(i,:)) ./ normalization_const;
                    allele_freq_power_grid_normalized = normalize_hist(x_grid, allele_freq_power_grid);
                    
                    figure; hold on; plot(x_grid, var_exp_grid, '.');
                    plot(x_grid, var_exp_power_grid, 'r.');
                    plot(x_grid, var_exp_power_grid_normalized, 'g.');
                    plot(x_vec, 2.*beta_vec.^2.*x_vec.*(1-x_vec), 'k*');
                    title('true and powered var. explained distribution');
                    xlabel('MAF'); ylabel('V'); legend('all var', 'observed var', 'observed var (normalized)');
                    figure; hold on; plot(x_grid, allele_freq_grid, '.');
                    plot(x_grid, allele_freq_power_grid, 'r.');
                    plot(x_grid, allele_freq_power_grid_normalized, 'g.');
                    plot(x_vec, 2.*beta_vec.^2.*x_vec.*(1-x_vec), 'k*');
                    title('true and powered allele freq distribution');
                    xlabel('MAF'); ylabel('V'); legend('all freq', 'observed freq', 'observed freq (normalized)');
                    ylim([0 20]);
                    
                    figure; plot(x_grid, exp(power_log_mat(i,:)), '.');
                    title('power for constant \beta as function of MAF');
                    xlabel('MAF'); ylabel('power');
                end
            end % if debug plots
            if(j == 1) % first interval is [0,1/2]
                %               switch correction_mode
                %                   case 'exact'
                pow_integral_ratio(i,ss_ind) = integral_hist(x_grid,  ...
                    var_explained_density_mat(ss_ind,:) .* exp(power_log_mat(i,:))) ./ ...
                    pow_integral_denominator(ss_ind); % take power over all the way [0,1/2]
                %                   case 'noisy' % here we deal also with variance in estimating beta
                pow_integral_ratio_noisy(i,ss_ind) = integral_hist(x_grid,  ...
                    2.*var_explained_density_mat(ss_ind,:) .* vec2row( ...
                    (sigma_beta(i)^2 + noisy_beta_vec(i)^2) .* ...
                    (1 - Phi((beta_discovery_grid-noisy_beta_vec(i))./sigma_beta(i))) + ...
                    sigma_beta(i) * (noisy_beta_vec(i) + beta_discovery_grid) .* ...
                    exp(-(noisy_beta_vec(i) - beta_discovery_grid).^2 ./ (2.*sigma_beta(i).^2)) ./ ...
                    (sqrt(2*pi)))) ./ ...
                    pow_integral_denominator(ss_ind); % take power over all the way [0,1/2]. This comes from the two-dimensional integral !!!! 
                %                end
            end
            if(b_ind > a_ind) % non empty interval for correction
                % New: V -> V *(1 + pi_not(a,b) / pi(0,1/2))
                %            pow_integral(i,j) = 1 +
                %                switch correction_mode
                %                    case 'exact'
                pow_integral_one_minus_mat(i,j,ss_ind) = integral_hist(x_grid(a_ind:b_ind),  ...
                    var_explained_density_mat(ss_ind,a_ind:b_ind) .* (1-exp(power_log_mat(i,a_ind:b_ind)))) ./ ...
                    pow_integral_denominator(ss_ind); % take one minus power. Denominator goes all the way [0,1/2]
                %                    case 'noisy' % here we deal also with variance in estimating beta
                pow_integral_one_minus_mat_noisy(i,j,ss_ind) = integral_hist(x_grid(a_ind:b_ind),  ...
                    2.*var_explained_density_mat(ss_ind,a_ind:b_ind) .* vec2row( ...
                    (sigma_beta(i)^2 + noisy_beta_vec(i)^2) .* ...
                    Phi((beta_discovery_grid(a_ind:b_ind)-noisy_beta_vec(i))./sigma_beta(i)) - ...
                    sigma_beta(i) * (noisy_beta_vec(i) + beta_discovery_grid(a_ind:b_ind)) .* ...
                    exp(-(noisy_beta_vec(i) - beta_discovery_grid(a_ind:b_ind)).^2 ./ (2.*sigma_beta(i).^2)) ./ ...
                    (sqrt(2*pi)))) ./ ...
                    pow_integral_denominator(ss_ind); % take one minus power. Denominator goes all the way [0,1/2]
                
                
                pow_integral_only(i,j) = integral_hist(x_grid(a_ind:b_ind), ...
                    exp(power_log_mat(i,a_ind:b_ind)));
                %                end
            end
        end % loop on s index values
    end % loop on beta values % loop on s_ind
end % loop on [a,b] limit interval

pow_integral_one_minus = zeros(num_loci, num_corrections);
pow_integral_one_minus_noisy = pow_integral_one_minus;
for j=1:num_corrections
    pow_integral_one_minus(:,j) = pow_integral_one_minus_mat(:,j,s_ind(j));
    pow_integral_one_minus_noisy(:,j) = pow_integral_one_minus_mat_noisy(:,j,s_ind(j));
end

V_corrected = repmat(V, 1, num_corrections) .* ...
    (1 + pow_integral_one_minus ./ ...
    repmat(pow_integral_ratio(:,s_ind(1)), 1, num_corrections)); % V * (1 + \bar{pi}(beta,s,a,b) / pi(beta,s))
V_corrected = min(V_corrected, repmat(V.*MAX_CORRECTION, 1, num_corrections)); % threshold loci with too much power correction 

%V_corrected = V_corrected(:,1:4); % not using the last one for now. Why???

V_corrected_noisy = repmat(V, 1, num_corrections) .* ...
    (1 + pow_integral_one_minus_noisy ./ ...
    repmat(pow_integral_ratio_noisy(:,s_ind(1)), 1, num_corrections)); % V * (1 + \bar{pi}(beta,s,a,b) / pi(beta,s))
V_corrected_noisy = V_corrected_noisy(:,1:4); % not using the last one for now

V_corrected_mat = zeros(num_corrections, length(s_grid), num_loci);
for j=1:num_corrections
    V_corrected_mat(j,:,:) = repmat(vec2row(V), length(s_grid), 1) .* ( 1 + ...
        reshape(pow_integral_one_minus_mat(:,j,:), length(beta_vec), length(s_grid))' ./ ...
        reshape(repmat(pow_integral_ratio, 1, 1), length(beta_vec), length(s_grid))');
end

% Naive is pointwise correction (doesn't include pop. gen.???)
V_corrected_naive = vec2column(V ./ vec2column(exp(power_log_vec))); % just divide by power (true)
V_corrected_strings_naive{1} = 'Park (true-effect)';
V_corrected_naive(:,2) = V_noisy ./ vec2column(exp(noisy_power_log_vec)); % just divide by power (observed)
V_corrected_strings_naive{2} = 'Park (observed-effect)';

% pointwise correction for integrating over beta
for j=1:num_loci
    V_corrected_naive(j,3) = V_noisy(j) ./ ... % vec2column(exp(noisy_power_log_vec));
        integral_hist(beta_grid, exp(power_log_beta_mat(j,:)) .* observed_beta_hist(j,:));
end % divide by average power (observed)
V_corrected_strings_naive{3} = 'Park (observed-effect-var-beta)';

V_corrected_naive(:,4) = V .* x_vec ./ pow_integral_only(:,1); % correction by Eliana (assuming s=0. Perform integral from zero to beta for each locus)
V_corrected_strings_naive{4} = 'Pop. Gen. [0,f] s=0';
V_corrected_naive(:,5) = V .* (x_vec-0.01) ./ pow_integral_only(:,5); % correction by Eliana (assuming s=0. Perform integral from zero to beta for each locus)
V_corrected_strings_naive{5} = 'Pop. Gen. [0.01,f] s=0'; % common loci
V_corrected_naive = min(V_corrected_naive, repmat(V.*2000000000, 1, size(V_corrected_naive,2))); % threshold loci with too much power correction 

V_corrected = [V V_corrected V_corrected_naive]; % add naive pointwise correction (this makes it not comparable to power !) 
V_corrected_strings = ['Observed' vec2row(V_corrected_strings) V_corrected_strings_naive];
