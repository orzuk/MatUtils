% A script for finding the maximum likelihood parameters given a set of found loci
%function [s a] = test_maximize_allele_freq_likelihood(x_vec, true_beta_vec, N, mu, ...
%    num_cases, num_controls, alpha)

%load(??); % load data
ttt = cputime;
AssignGeneralConstants;
new_loci = 1; % generate new loci or just re-analyze
test_beta_observation = 0; % test empirical beta estimation
plot_moments=1; %  just plot theoretical moments
run_simulations = 0;  % run simulations (heavy part)
if(~exist('html_outdir', 'var'))
    [machine machine_delim html_outdir] = get_machine_type();
end
test_type = 'single-locus'; test_stat = 'chi-square-QTL'; % default (QTLs for now)

disease_html_root_dir = fullfile(html_outdir, 'common_disease_broad_catalog');
simulate_loci = 0; % use real or simulated data
use_power = 1; % whether to include power calculation in model

constant_bin_number = 0; % 0 - equal spacing of bins, 1 equal # of points in a bin
figs_dir = '../../common_disease_model/figs';
%num_corrections = 5;
N=10^4; x_grid = 0.5.*(1:N)./N; % set effective population size and MAF grid
alpha = 5*10^(-8); % set significance threshold


% load(fullfile(disease_html_root_dir, 'disease_data.mat')); % load data
load('disease_data.mat'); num_snps = length(data.GRR);
H = load('../../common_disease_model/data/Height/all_height_loci_effect_sizes.mat');
height_inds = strmatch('Height', data.Trait);
[height_SNPS new_inds height_inds2] = intersect(data.SNPs(height_inds), H.SNP);
height_inds = height_inds(new_inds);
%data.SNPs(height_inds) H.SNP(height_inds2)
data.discovery_beta = zeros(num_snps,1);
data.discovery_beta(height_inds) = H.Beta(height_inds2);
data.combined_beta(height_inds) = H.BetaCombined(height_inds2);
find(data.Beta(height_inds) - H.BetaReplication(height_inds2))

special_traits = {'Height', 'Lipid HDL', 'Lipid LDL', 'Lipid TriGly', ...
    'Breast Cancer', 'Crohn''s disease', 'Type 1 diabetes', 'Type 2 diabetes'};
alpha_vec = [5*10^(-8) 5*10^(-8) 5*10^(-8) 5*10^(-8) ...
    5*10^(-8) 5*10^(-8) 5*10^(-8) 5*10^(-8)]; % first one is 10^(-8) or 10^(-6) ??
fit_mode = 'selection'; simulate_mode = 'slope'; % 'selection'; % 'slope'; % 'selection'; % 'slope';
my_fig_file = '../fig/simulated_s_0';

if(run_simulations)
    for t=[1] %   [5 7 8] % 1:1 % length(special_traits) % loop on traits. Now just do Height (1) and Crohn's (6)
        alpha = alpha_vec(t);
        if(~simulate_loci)
            trait_is = strrep(special_traits{t}, ' ', '-')
        else
            trait_is = ['simulated'];
        end
        disease_ind = strmatch(special_traits{t}, data_params.trait_name);
        snp_inds = strmatch(special_traits{t}, data.Trait);
        if(~isempty(snp_inds))
            if(new_loci)
                mu = 0.05;
                num_cases =  [data.discovery_num_cases(snp_inds(1)) ...
                    data.replication_num_cases(snp_inds(1))];
                num_controls = [data.discovery_num_controls(snp_inds(1)) ...
                    data.replication_num_controls(snp_inds(1))];
                if(max(num_controls) < 0)
                    num_controls = [];
                end
                switch data.Trait{snp_inds(1)} % special correction: replication add sample size
                    case 'Height' % use replication set size. Why???
                        %                        num_cases = 183727; % 133653; % 133653 is stage 1 sample size % 183727 is combined sample size
                    case 'Crohn''s disease'
                        num_cases = [6333 15694]; % set discovery and replication stages
                        num_controls = [15056 14026];
                        alpha = [5*10^(-6) 0.01]; % set significance thresholds for discovery and replication stages
                end
                
                num_loci = length(snp_inds);
                theta = 0; s = [];
                s_sim = -0.0001; a_sim = -0.001;
                switch simulate_mode
                    case 'selection'
                        s_or_a_sim = s_sim;
                    case 'slope'
                        s_or_a_sim = a_sim;
                end
                % num_cases = 20000;
                
                if(simulate_loci)
                    num_cases = 10000; % change to adjust average power
                    num_loci = 50; % increase data size (only a few pass discovery)
                    [ true_beta_vec x_vec grr_vec h_liab_vec ...
                        discovery_indicator_vec discovery_power_vec ...
                        noisy_beta_vec replication_beta_vec] = ...
                        simulate_loci_effect_size(1, N, mu, theta, ...
                        simulate_mode, s_or_a_sim, num_cases, num_controls, alpha, ...
                        num_loci, use_power);
                    discovery_effect_flag = 1; % we always have all effect sizes
                    trait_type = 'QTL'; % for now simualte just QTLs
                else
                    x_vec = data.MAF(snp_inds);
                    true_beta_vec = data.Beta(snp_inds);
                    true_grr_vec = data.GRR(snp_inds);
                    true_mu = data.Prevalence(snp_inds(1)); mu = true_mu;
                    trait_type = data.trait_type{snp_inds(1)}
                    
                    
                    switch lower(data.trait_type{snp_inds(1)}) % convension: always take MAF !
                        case 'binary'
                            [x_vec true_grr_vec] = flip_allele(x_vec, true_grr_vec, 'binary', 1);
                            true_beta_vec = [true_grr_vec repmat(true_mu, num_loci, 1)];
                            test_type = 'armitage';
                            test_stat = 'chi-square';
                        case 'qtl'
                            [x_vec true_beta_vec] = flip_allele(x_vec, true_beta_vec, 'QTL', 1);
                            test_type = 'single-locus';
                            test_stat = 'chi-square-QTL';
                            num_controls = [];
                    end
                    if( (isfield(data, 'discovery_beta')) && (any(data.discovery_beta(snp_inds))) )
                        noisy_beta_vec = vec2column(data.discovery_beta(snp_inds));
                        discovery_effect_flag = 1; % we have discovery effect sizes
                    else
                        noisy_beta_vec = true_beta_vec; % we actually don't know the true ones ...
                        discovery_effect_flag = 0; % we don't have discovery effect sizes
                    end
                    replication_beta_vec = true_beta_vec;
                    discovery_indicator_vec = ones(size(x_vec,1),1);
                    discovery_power_vec = ones(size(x_vec,1),1);
                end
                
                switch length(num_cases)
                    case 1 % just one study size. Use it for all
                        num_cases_discovery = num_cases;
                        num_cases_replication = num_cases;
                        num_cases_combined = num_cases;
                        num_controls_discovery = num_controls;
                        num_controls_replication = num_controls;
                        num_controls_combined = num_controls;
                    case 2 % discovery and replication
                        num_cases_discovery = num_cases(1);
                        num_cases_replication = num_cases(2);
                        num_cases_combined = sum(num_cases);
                        switch trait_type
                            case {'binary', 'Binary'}
                                num_controls_discovery = num_controls(1);
                                num_controls_replication = num_controls(2);
                                num_controls_combined = sum(num_controls);
                            case {'QTL', 'Quantitative'}
                                num_controls_discovery = [];
                                num_controls_replication = [];
                                num_controls_combined = [];
                        end
                end
                
                [beta_discovery_grid beta_discovery_inv_grid] = ...
                    compute_discovery_boundary(x_grid, mu, ...
                    num_cases_discovery, num_controls_discovery, alpha(1), ...
                    test_type, test_stat); % Determine discovery boundry
                
                true_beta_vec = abs(true_beta_vec); % make everything positive (easier)
                num_observed_loci = sum(discovery_indicator_vec);
                noisy_beta_vec = abs(noisy_beta_vec); % make sure they're all positive
                [noisy_beta_vec sort_perm] = sort(noisy_beta_vec(:,1));
                true_beta_vec = true_beta_vec(sort_perm,:);
                replication_beta_vec = replication_beta_vec(sort_perm,:);
                x_vec = x_vec(sort_perm); % sort by effect sizes
                discovery_indicator_vec = discovery_indicator_vec(sort_perm); discovery_inds_vec = find(discovery_indicator_vec);
                discovery_power_vec = discovery_power_vec(sort_perm);
                x_MAF_vec = min(x_vec, 1-x_vec);
                x_observed_vec = x_vec(discovery_inds_vec);
                true_beta_observed_vec = true_beta_vec(discovery_inds_vec,:);
                noisy_beta_observed_vec = noisy_beta_vec(discovery_inds_vec,:);
                replication_beta_observed_vec = replication_beta_vec(discovery_inds_vec,:);
                discovery_power_observed_vec = discovery_power_vec(discovery_inds_vec);
                
                if(constant_bin_number)
                    bin_size = 20; % num_loci; % 18; % num_loci; % 18; % num_loci; % length(x_vec) / num_bins;
                    %            num_bins = floor(length(x_vec) / bin_size); % 18; % Divide to bins
                    [bin_observed_starts bin_observed_ends bin_lengths num_bins] = ... % bin_ind] = ...
                        divide_region_to_blocks(1, length(x_observed_vec), bin_size); % blocks of equal size
                else % here constant bin size
                    num_bins = 1;
                    [hhh, bin_boundaries] = hist(noisy_beta_vec(:,1), num_bins);
                    if(length(bin_boundaries) > 1)
                        bin_size = bin_boundaries(2)-bin_boundaries(1);
                    else
                        bin_size = noisy_beta_vec(end) - noisy_beta_vec(1) + epsilon;
                    end
                    bin_boundaries = ...
                        [min(noisy_beta_vec(1,1), min(true_beta_observed_vec(:,1)))-epsilon ...
                        bin_boundaries+0.5*bin_size+epsilon];
                    bin_boundaries(end) = max(bin_boundaries(end), ...
                        max(true_beta_observed_vec(:,1) + epsilon)); % make sure true effect sizes are within range
                    [~, all_bin_inds] = histc(true_beta_observed_vec(:,1), bin_boundaries);
                    observed_bins = unique(all_bin_inds); % bins with observed points
                    empty_bins = setdiff(1:num_bins, observed_bins);
                    observed_bins_indicator = zeros(num_bins,1);
                    observed_bins_indicator(observed_bins) = 1;
                    bin_observed_starts = zeros(num_bins, 1); bin_observed_ends = zeros(num_bins, 1);
                    bin_observed_starts(observed_bins) = [1 vec2row(find(diff(all_bin_inds))+1)];
                    bin_observed_ends(observed_bins) = [vec2row(find(diff(all_bin_inds))) num_observed_loci];
                    
                    [~, all_bin_inds] = histc(noisy_beta_vec(:,1), bin_boundaries);
                    occupied_bins = unique(all_bin_inds); % bins with observed points
                    bin_starts = zeros(num_bins, 1); bin_ends = zeros(num_bins, 1);
                    bin_starts(occupied_bins) = [1 vec2row(find(diff(all_bin_inds))+1)];
                    bin_ends(occupied_bins) = [vec2row(find(diff(all_bin_inds))) num_loci];
                end
                
                figure; hold on; legend_vec = [];
                if(simulate_loci) % here we know the true values
                    plot(x_MAF_vec, abs(true_beta_vec), '.');
                    plot(x_MAF_vec(discovery_inds_vec), ...
                        abs(true_beta_vec(discovery_inds_vec)), 'rx');
                    legend_vec = [legend_vec {'missed (true)', 'discovered (true)'}];
                end
                if((exist('noisy_beta_vec', 'var')) && discovery_effect_flag) % this is the DISCOVERY effect size
                    if(simulate_loci) % here we know what we've missed
                        plot(x_MAF_vec, abs(noisy_beta_vec), 'o');
                        legend_vec = [legend_vec {'missed (discovery)'}];
                    end
                    plot(x_MAF_vec(discovery_inds_vec), ...
                        abs(noisy_beta_vec(discovery_inds_vec)), 'ro');
                    legend_vec = [legend_vec {'discovered (discovery)'}];
                end
                if(exist('replication_beta_vec', 'var')) % this is the REPLICATION effect size
                    if(simulate_loci) % here we know what we've missed
                        plot(x_MAF_vec, abs(replication_beta_vec), '+');
                        legend_vec = [legend_vec {'missed (replication)'}];
                    end
                    plot(x_MAF_vec(discovery_inds_vec), ...
                        abs(replication_beta_vec(discovery_inds_vec,1)), 'r+');
                    legend_vec = [legend_vec {'discovered (replication)'}];
                end
                
                plot(x_grid, [beta_discovery_grid beta_discovery_inv_grid], 'g--', 'linewidth', 3); % plot discovery boundary
                title(['Effect Size vs. Allele. Freq. All Loci, ' trait_is]);
                xlabel('MAF');
                legend_vec = [legend_vec {'discovery-boundary'}];
                legend(legend_vec);
                for j=1:length(bin_boundaries)
                    line([0 0.5], [bin_boundaries(j) bin_boundaries(j)], ...
                        'color', 'k', 'linestyle', '--');
                end
                switch trait_type
                    case {'Binary', 'binary'}
                        plot(x_grid, repmat(1, length(x_grid), 1), 'k:'); % plot line GRR=1
                        y_min = 0; effect_str = 'GRR';
                    case {'QTL', 'Quantitative'}
                        y_min = 0; effect_str = '\beta';
                end
                ylabel(effect_str);
                ylim([y_min min(max(max(abs(noisy_beta_vec(:,1)), abs(true_beta_vec(:,1))))*1.1, 50)]); % make scale smaller
                my_saveas(gcf, fullfile(figs_dir, trait_is, ['beta_vs_MAF_' trait_is]), format_fig_vec);
            end % if new loci
            for s_input_ind = 1:1 % try fitted (2) or true (1) s for variance correction
                if(new_loci)
                    if(s_input_ind == 1)
                        s_input = s_sim; s_str = ' true s';
                    else
                        s_input = []; s_str = ' fitted s';
                    end
                    mean_beta = zeros(1, num_bins);
                    V_explained = zeros(1, num_bins);
                    V_true = zeros(1, num_bins);
                    V_true_common = zeros(1, num_bins);
                    V_true_rare = zeros(1, num_bins);
                    first_bin_flag =  1;
                    s = cell(num_bins,1);
                    a = cell(num_bins,1);
                    LL = cell(num_bins,1);
                    V = cell(num_bins,1);
                    V_corrected = cell(num_bins,1);
                    V_corrected_noisy = cell(num_bins,1);
                    V_true_vec = cell(num_bins,1);
                    pow_integral = cell(num_bins,1);
                    power_log_vec = cell(num_bins,1);
                    V_corrected_strings = cell(num_bins,1);
                    for i=1:num_bins % loop on effects of different sizes
                        run_bin = i
                        bin_inds = bin_starts(i):bin_ends(i); % ((i-1)*bin_size+1):i*bin_size;
                        if(~ismember(0, bin_inds))
                            V_true_vec{i} = 2 .* x_vec(bin_inds) .* (1-x_vec(bin_inds)) .* ...
                                true_beta_vec(bin_inds,1).^2;
                            V_true(i) = sum(V_true_vec{i});
                            common_inds = find(x_MAF_vec(bin_inds) > 0.01);
                            rare_inds = find(x_MAF_vec(bin_inds) <= 0.01);
                            V_true_common(i) = sum(V_true_vec{i}(common_inds));
                            V_true_rare(i) = sum(V_true_vec{i}(rare_inds));
                        end
                        if(observed_bins_indicator(i))
                            bin_observed_inds = bin_observed_starts(i):bin_observed_ends(i); % ((i-1)*bin_size+1):i*bin_size;
                            [s{i} a{i} LL{i} V{i} V_corrected{i} V_corrected_noisy{i} ...
                                V_corrected_strings{i} pow_integral{i} power_log_vec{i}] = ...
                                maximize_gwas_observed_allele_freq_and_effect_size_likelihood( ...
                                x_observed_vec(bin_observed_inds), ...
                                true_beta_observed_vec(bin_observed_inds,:), ...
                                replication_beta_observed_vec(bin_observed_inds,:),...
                                x_grid, beta_discovery_grid, ...
                                N, mu, num_cases, num_controls, alpha, use_power, ...
                                s_input, 'replication', trait_type, 1); % find optimal s, (beta) and compute correction
                            if(first_bin_flag)
                                num_corrections =  size(V_corrected{i},2);
                                V_explained_corrected = zeros(num_bins, num_corrections); % currently use 4 differnet corrections
                                first_bin_flag = 0;
                                if(simulate_loci)
                                    my_fig_file = ['../fig/simulated_s_' num2str(s_sim)];
                                    V_corrected_strings{i} = [V_corrected_strings{i} ...
                                        'true', 'true-common>0.01'];
                                    show_corrections = {'Observed', ...
                                        'Pop. Gen. [0, 0.5], fitted s', ...
                                        'Pop. Gen. [0, 0.5], s=0', ...
                                        'Park (observed-effect)', ...
                                        'Park (observed-effect-var-beta)', ...
                                        'Pop. Gen. [0, f], s=0', ...
                                        'true', 'true-common>0.01'};
                                    
                                    % % %                                     'Pop. Gen. [0, f], fitted s', ...
                                    % % %                                         'Pop. Gen. [0.01, 0.5], fitted s', ...
                                    % % %                                         'Park (true-effect)', ...
                                    % % %                                         'Pop. Gen. [0.01, f], s=0', ...
                                    % % %                                         'Pop. Gen. [0.01, f], fitted s', ...
                                    
                                    
                                    %                                show_inds = [1 4 3 6 7 num_corrections+1 num_corrections+2]; % what to plot in variance explained
                                else
                                    my_fig_file = ['../fig/' special_traits{t}];
                                    show_corrections = {'Observed', ...
                                        'Pop. Gen. [0, f], fitted s', ...
                                        'Pop. Gen. [0, 0.5], fitted s', ...
                                        'Pop. Gen. [0, 0.5], s=0', ...
                                        'Park (observed-effect)', ...
                                        'Park (observed-effect-var-beta)', ...
                                        'Pop. Gen. [0, f], s=0'};
                                    %
                                    %                                                                             'Pop. Gen. [0.01, 0.5]', ...
                                    %                                         'Pop. Gen. [0.01, f]', ...
                                    %                                         'Pop. Gen. [0.01, f] s=0'};
                                    
                                    %                                show_inds = [1 4 3 6 7];  % what to plot in variance explained (no true values)
                                end
                                [~, show_inds] = intersect(V_corrected_strings{1}, ...
                                    show_corrections);
                                
                            end % if first bin
                            
                            discovery_power_is = discovery_power_observed_vec(bin_observed_inds)
                            my_saveas(gcf, fullfile(figs_dir, trait_is, ...
                                ['var_explained_vs_s_sensitivity']), format_fig_vec);
                            
                            %                        if(num_bins == 1)
                            %                           my_saveas(gcf, my_fig_file, format_fig_vec);
                            %                        end
                            mean_beta(i) = mean(true_beta_vec(bin_inds));
                            V_explained(i) = sum(V{i});
                            V_explained_corrected(i,:) = sum(V_corrected{i},1); % get orientation right
                            
                            %                if(i == 1) % prepare a table
                            tmp_pow = zeros(length(bin_inds),1);
                            tmp_pow(find(discovery_indicator_vec(bin_inds))) = pow_integral{i}(:,1);
                            R = [vec2column(true_beta_vec(bin_inds)) vec2column(x_MAF_vec(bin_inds)) ...
                                vec2column(discovery_power_vec(bin_inds)) vec2column(discovery_indicator_vec(bin_inds)) ...
                                vec2column(V_true_vec{i}) vec2column(V_true_vec{i}) .* vec2column(discovery_indicator_vec(bin_inds)) tmp_pow];
                            R = [R' sum(R,1)']'; R(end,end) = R(end,end-1) / R(end,end-2); % sum and also get true ratio
                            column_labels = {'beta', 'MAF', 'Power', 'Discovered', 'V^2', 'V^2 discovered', 'Power-var-integral'};
                            R = [column_labels' num2cell(R)']'; R = num2str_cell(R, 5);
                            savecellfile(R, ['debug_power_correction_bin' num2str(i) '.txt']);
                            %                end
                            
                            figure; hold on; plot(x_MAF_vec(bin_inds), ...
                                abs(true_beta_vec(bin_inds)), '.'); % simple plot of all loci
                            plot(x_MAF_vec(bin_inds(find(discovery_indicator_vec(bin_inds)))), ...
                                abs(true_beta_vec(bin_inds(find(discovery_indicator_vec(bin_inds))))), 'r*');
                            title(['Effect Size vs. Allele. Freq. Bin #' num2str(i)]);
                            xlabel('MAF'); ylabel(effect_str);
                            legend({'missed', 'discovered'});
                        end
                    end % loop on bins
                end % if new loci
                
                switch trait_type
                    case {'QTL', 'Quantitative'}
                        beta_grid = -1:0.001:1; % all possible betas
                    case 'Binary'
                        beta_grid = 0:0.001:10; % all possible GRR's
                end
                
                switch fit_mode % many different plots
                    case 'slope'
                        figure; hold on; plot(mean_beta, a, '.');
                        title('Selection Slope vs. Effect Size');
                        xlabel(effect_str); ylabel('slope');
                        if(simulate_loci) % compare to true values
                            plot(mean_beta, a_sim .* mean_beta, 'r');
                            legend('fitted', 'true');
                        end
                        
                    case 'selection'
                        figure; hold on; plot(mean_beta(observed_bins), ...
                            cell2vec(s(observed_bins)), '.');
                        title('Selection Coefficient vs. Effect Size');
                        xlabel(effect_str); ylabel('s');
                        lin_fit = polyfit(mean_beta(observed_bins), ...
                            cell2vec(s(observed_bins)), 1);
                        plot(mean_beta, lin_fit(1) .* mean_beta + lin_fit(2), 'k')
                        if(simulate_loci) % compare to true values
                            switch simulate_mode
                                case 'selection'
                                    plot(mean_beta, repmat(s_sim, length(mean_beta), 1), 'r');
                                case 'slope'
                                    plot(mean_beta, a_sim .* mean_beta, 'r');
                            end
                            legend('fitted-s', 'linear-fit', 'true');
                        else
                            legend('fitted-s', 'linear-fit');
                        end
                        my_saveas(gcf, fullfile(figs_dir, trait_is, ['beta_vs_s']), format_fig_vec);
                        
                        
                        %                    figure; hold on; % Bar plot of var. explained vs. effect size
                        %                    bar(mean_beta, V_explained);
                        if(exist('bin_boundaries', 'var'))
                            [h bins] = weighted_hist(vec2row(true_beta_observed_vec(:,1)), ...
                                vec2row(cell2vec(V',[],0)), bin_boundaries);
                            
                        else
                            [h bins] = weighted_hist(vec2row(true_beta_observed_vec(:,1)), ...
                                vec2row(cell2vec(V',[],0)), num_bins-1);
                        end
                        V_corrected_vec = cell2vec(V_corrected(observed_bins),[],0);
                        corrected_h = zeros(size(V_explained_corrected, 2), num_bins);
                        for j=1:size(V_explained_corrected, 2) % loop on different corrections
                            if(exist('bin_boundaries', 'var'))
                                [corrected_h(j,:) corrected_bins{j}] = ...
                                    weighted_hist(vec2row(true_beta_observed_vec(:,1)), ...
                                    vec2row(V_corrected_vec(:,j)), bin_boundaries);
                            else
                                [corrected_h(j,:) corrected_bins{j}] = ...
                                    weighted_hist(vec2row(true_beta_observed_vec(:,1)), ...
                                    vec2row(V_corrected_vec(:,j)), num_bins-1);
                            end
                        end
                        num_bars = 2 + num_corrections;
                        corrected_h(j+1,:) = vec2row(V_true); % add true variance
                        corrected_h(j+2,:) = vec2row(V_true_common); % add true common variance
                        
                        total_V_explained_integral = zeros(length(show_inds),1);
                        for jj=1:length(show_inds) % num_corrections. Compute var. explaine for each one
                            if(show_inds(jj) <= num_corrections)
                                total_V_explained_integral(jj) = corrected_h(show_inds(jj));
                                %                                ...
                                %                                    integral_hist(beta_grid, ...
                                %                                    all_beta_estimate_corrected(show_inds(jj),:));
                            else % here we deal with true vec (???)
                                bin_locs = ...
                                    bin_boundaries(1):range(bin_boundaries)/1000 ...
                                    :bin_boundaries(end);
                                if(~isempty(strfind(V_corrected_strings{1}{show_inds(jj)}, ...
                                        'common')))
                                    include_var_inds = common_inds;
                                else
                                    include_var_inds = 1:length(true_beta_vec);
                                end
                                total_V_explained_integral(jj) = sum(V_true_vec{1}(include_var_inds));    %  total_V_explained_integral(1) = integral_hist(beta_grid, all_beta_estimate);
                            end
                        end % loop on difference corrections
                        total_V_explained_str = ... % show var. explained in %
                            num2str_cell(num2cell(total_V_explained_integral*100),3,[],1);
                        legend_vec = V_corrected_strings{1}; %, ...
                        % 'true' 'true-common>0.01']; % new legends
                        %                        legend_vec = {'observed', 'corrected<f^*', ...
                        %                            'corrected>0.01', 'corrected', ...
                        %                            'corrected(neutral,s=0)', 'corrected-pointwise-true', 'corrected-pointwise-noisy', ...
                        %                            '', 'true', 'true-common>0.01'};
                        cur_legend_vec = cell(length(show_inds),1);
                        for jj = 1:length(show_inds)
                            cur_legend_vec{jj} = [legend_vec{show_inds(jj)} ...
                                ', h^2=' total_V_explained_str{jj}];
                        end
                        figure; hold on; % Figure: bars ...
                        if(num_bins > 1)
                            show_bar = vec2row(corrected_h);
                            bar(bins, show_bar(:,show_inds)); % bar_multi
                        else
                            show_bar = [zeros(num_bars, 1) vec2column(corrected_h)]';
                            bar([bins-epsilon bins], show_bar(:,show_inds)); % bar_multi
                        end
                        title(['Var. Explained for different ' effect_str ' ' s_str]);
                        xlabel(effect_str); ylabel('Var. Explained');
                        legend(cur_legend_vec);
                        my_saveas(gcf, fullfile(figs_dir, trait_is, ...
                            'var_explained_by_bin'), format_fig_vec);
                        figure; hold on; % figure ... bars ...
                        explained_sum = sum(corrected_h',1);
                        show_bar = [zeros(num_bars, 1) explained_sum']';
                        bar(show_bar(:,show_inds));  % bar sum
                        title(['Total Var. Explained for all different ' effect_str ' ' s_str]);
                        ylabel('Var. Explained');
                        legend(cur_legend_vec);
                        my_saveas(gcf, fullfile(figs_dir, trait_is, 'var_explained_total'), format_fig_vec);
                        
                        
                        figure; hold on; % New: figure of var. explained dependency on beta
                        %                        all_beta_estimate = zeros(size(beta_grid));
                        all_beta_estimate_corrected = zeros(num_corrections, length(beta_grid));
                        for j=1:num_observed_loci % loop and add contribution to variance explained
                            cur_replication_beta = replication_beta_observed_vec(j,:);
                            cur_x = x_observed_vec(j);
                            
                            switch trait_type
                                case {'QTL', 'Quantitative'}
                                case {'Binary', 'binary'}
                                    if(cur_replication_beta(1) < 1)
                                        cur_replication_beta(1) = 1 / cur_replication_beta(1);
                                        cur_x = 1-cur_x;
                                    end
                            end
                            cur_true_beta_estimate = ... % estimation of true beta (indep. of correction used)
                                observed_effect_size_dist(cur_replication_beta, ...
                                cur_x, test_stat, alpha, beta_grid, ...
                                num_cases_replication, num_controls_replication, 'replication');
                            %                            all_beta_estimate = all_beta_estimate + ...
                            %                                cur_true_beta_estimate .* V{1}(j); % observed variance
                            for jj=1:num_corrections % look how much is corrected for each correction type
                                all_beta_estimate_corrected(jj,:) = ...
                                    all_beta_estimate_corrected(jj,:) + ...
                                    cur_true_beta_estimate .* V_corrected{1}(j,jj); % corrected variance (pointwise)
                            end
                        end
                        switch trait_type % collapse MAF [0,0.5] and [0.5,1]
                            case 'QTL' % flip beta to collapse positive and negative together
                                mid_ind = (length(beta_grid)-1)/2;
                                beta_grid = beta_grid(mid_ind+2:end);
                                all_beta_estimate_corrected = all_beta_estimate_corrected(:,mid_ind+2:end) + ...
                                    all_beta_estimate_corrected(:,mid_ind:-1:1);
                            case {'Binary', 'binary'}
                                mid_ind = find(beta_grid < 1, 1, 'last');
                                %                                beta_inv_grid = min(100, 1 ./ beta_grid(1:mid_ind)); % don't get too high GRR
                                beta_grid = beta_grid(mid_ind+2:end);
                                all_beta_estimate_corrected = ...
                                    all_beta_estimate_corrected(:,mid_ind+2:end);
                                %                                for jj=1:num_corrections
                                %                                    [temp_beta_grid temp_all_beta_estimate_corrected(jj,:)] = ...
                                %                                    sum_hist(beta_grid, all_beta_estimate_corrected(jj,mid_ind+2:end), ...
                                %                                    beta_inv_grid(1:mid_ind), ...
                                %                                    all_beta_estimate_corrected(jj,1:mid_ind), [], 0);
                                %                                end
                                %                                beta_grid = temp_beta_grid;
                                %                                all_beta_estimate_corrected = temp_all_beta_estimate_corrected;
                        end
                        
                        %                        plot(beta_grid, all_beta_estimate);
                        line_width_vec = ones(length(show_inds), 1); line_width_vec(1) = 3; % make true effect size thicker
                        for jj=1:length(show_inds) % num_corrections
                            if(show_inds(jj) <= num_corrections)
                                plot(beta_grid, all_beta_estimate_corrected(show_inds(jj),:), ...
                                    color_vec(jj), 'linewidth', line_width_vec(jj));
                                %                              total_V_explained_integral(jj) = ...
                                %                                  integral_hist(beta_grid, ...
                                %                                  all_beta_estimate_corrected(show_inds(jj),:));
                            else % here we deal with true vec (???)
                                bin_locs = ...
                                    bin_boundaries(1):range(bin_boundaries)/1000 ...
                                    :bin_boundaries(end);
                                if(~isempty(strfind(V_corrected_strings{1}{show_inds(jj)}, ...
                                        'common')))
                                    include_var_inds = common_inds;
                                else
                                    include_var_inds = 1:length(true_beta_vec);
                                end
                                %                               total_V_explained_integral(jj) = sum(V_true_vec{1}(include_var_inds));    %  total_V_explained_integral(1) = integral_hist(beta_grid, all_beta_estimate);
                                
                                [hhh bbbins] = weighted_hist(vec2row(true_beta_vec(include_var_inds)), ...
                                    vec2row(V_true_vec{1}(include_var_inds)), bin_locs);
                                plot(bbbins(find(hhh)), hhh(find(hhh)), [color_vec(jj), 'x']);
                                %                                hist_density(true_beta_vec, V_true_vec{1}, 'k');
                            end
                        end % loop on difference corrections
                        title('Var. explained as function of effect size');
                        xlabel(effect_str); ylabel('Var. Explained'); % should be '|\beta|'
                        
                        legend(cur_legend_vec);
                        switch trait_type
                            case {'QTL', 'Quantitative'}
                                xlim([0 0.3]); % temp. effect size range ..
                            case {'Binary', 'binary'}
                                xlim([1 max(noisy_beta_vec)+1]);
                        end
                        my_saveas(gcf, fullfile(figs_dir, trait_is, ...
                            [trait_is '_var_density']), format_fig_vec);
                        
                        figure; hold on; % plot cummulative var explained as dependency on beta
                        %                        plot(beta_grid, cumsum(all_beta_estimate));
                        for jj=1:length(show_inds) % num_corrections
                            if(show_inds(jj) <= num_corrections)
                                plot(beta_grid, (beta_grid(2)-beta_grid(1)) .* ...
                                    cumsum(all_beta_estimate_corrected(show_inds(jj),:)), ...
                                    color_vec(jj), 'linewidth', line_width_vec(jj));
                            else % last two are true
                                if(~isempty(strfind(V_corrected_strings{1}{show_inds(jj)}, ...
                                        'common')))
                                    include_var_inds = common_inds;
                                else
                                    include_var_inds = 1:length(true_beta_vec);
                                end
                                plot([0 vec2row(true_beta_vec(include_var_inds)) max(beta_grid)] , ...
                                    [0 vec2row(cumsum(V_true_vec{1}(include_var_inds))) sum(V_true_vec{1}(include_var_inds))], ...
                                    [color_vec(jj) '--']);
                                
                            end
                        end
                        title('Cum. Var. explained as function of effect size');
                        xlabel(effect_str); ylabel('Var. Explained'); % should be '|\beta|'
                        legend(cur_legend_vec);
                        %                        legend(['observed ' total_V_explained_str{1}], ...
                        %                            ['corrected ' total_V_explained_str{2}]);
                        switch trait_type
                            case 'QTL'
                                xlim([0 0.3]);
                            case {'binary', 'Binary'}
                                xlim([1 max(noisy_beta_vec)+1]);
                        end
                        my_saveas(gcf, fullfile(figs_dir, trait_is, ...
                            [trait_is '_var_cumulative']), format_fig_vec);
                        
                end
            end % loop on s fitted vs. true
        end % isempty snp inds
    end % loop on different traits
end % if plot_moments


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s


if(test_beta_observation)
    true_beta = 0.12;
    f_vec = 0.05;
    trait_type = 'QTL';
    alpha = 5*10^(-4);
    beta_grid = -1:0.001:1;
    num_cases = 5000; % [10000 3000] ;
    num_controls = 12000; true_grr = 1.5; mu = 0.05; % set binary trait's parameters
    grr_grid = 0:0.001:5;
    study_type = 'replication'; % 'combined'; % 'discovery'; % 'replication';
    
    observed_beta_hist_sim = observed_effect_size_dist( true_beta, f_vec, trait_type, ...
        alpha, beta_grid, num_cases, num_controls, study_type, 'simulate');
    observed_beta_hist_analytic = observed_effect_size_dist( true_beta, f_vec, trait_type, ...
        alpha, beta_grid, num_cases, num_controls, study_type, 'gaussian');
    
    observed_grr_hist_sim = observed_effect_size_dist( [true_grr mu], f_vec, 'binary', ...
        alpha, grr_grid, num_cases, num_controls, study_type, 'simulate');
    observed_grr_hist_analytic = observed_effect_size_dist( [true_grr mu], f_vec, 'binary', ...
        alpha, grr_grid, num_cases, num_controls, study_type, 'gaussian');
    
    
    figure; hold on;
    x_lim = find( max(observed_beta_hist_sim, observed_beta_hist_analytic) > epsilon);
    plot(beta_grid, observed_beta_hist_sim, 'b');
    plot(beta_grid, observed_beta_hist_analytic, 'r');
    plot(true_beta, 0, 'g*');
    xlim([min(true_beta, beta_grid(x_lim(1))) beta_grid(x_lim(end))]);
    title('observed effect size distribution QTL');
    legend('Simulation', 'Gaussian approximation', ['true ' effect_str]);
    xlabel(effect_str); ylabel('freq.');
    
    figure; hold on;
    x_lim = find( max(observed_grr_hist_sim, observed_grr_hist_analytic) > epsilon);
    plot(grr_grid, observed_grr_hist_sim, 'b');
    plot(grr_grid, observed_grr_hist_analytic, 'r');
    plot(true_grr, 0, 'g*');
    xlim([min(true_grr, grr_grid(x_lim(1))) grr_grid(x_lim(end))]);
    title('observed effect size distribution (disease)');
    legend('Simulation', 'Gaussian approximation', 'true GRR');
    xlabel('GRR'); ylabel('freq.');
    
    
    %     figure; hold on; % plot also GRR on a log scale
    %     x_lim = find( max(observed_grr_hist_sim, observed_grr_hist_analytic) > epsilon);
    %     plot(log(grr_grid), observed_grr_hist_sim, 'b');
    %     plot(log(grr_grid), observed_grr_hist_analytic, 'r');
    %     plot(true_grr, 0, 'g*');
    %     xlim([min(true_grr, grr_grid(x_lim(1))) grr_grid(x_lim(end))]);
    %     title('observed effect size distribution (disease)');
    %     legend('Simulation', 'Gaussian approximation', 'true GRR');
    %     xlabel('GRR'); ylabel('freq.');
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
if(plot_moments) % Just plot moments of population and variance as function of beta (theoretical plot)
    quant_vec = [0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99];
    plot_var_explained_by_allele_freq_quantiles(quant_vec, N, figs_dir);
end % if plot moments

time_elapsed = cputime - ttt
