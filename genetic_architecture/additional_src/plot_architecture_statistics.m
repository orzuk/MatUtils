% Make various plots for all architectures
%
% Input:
% good_architectures - set of candidate architectures
% params_struct - parameters for each architecture
% good_architectures_plot_file - file where to save plots
%
% Output:
%
function plot_architecture_statistics(good_architectures, ...
    params_struct, good_architectures_plot_file)

seperate_plots = 1; % Whether to plot everything together or seperately
AssignGeneralConstants;
font_size = 9; % set the same font size for all figs (default is 9)
legend_font_size = 7;
precision = 3;

set(0,'DefaultAxesFontWeight','bold');
% format_fig_vec = {'fig', 'jpg'}; % no eps
if(~exist('good_architectures_plot_file', 'var') || isempty(good_architectures_plot_file))
    save_flag = 0;
else
    save_flag = 1;
end
if(save_flag)
    my_mkdir(dir_from_file_name(good_architectures_plot_file));
end

num_architectures = length(good_architectures); % find strings
if(num_architectures > 0) % found architectures
    N_vec = zeros(num_architectures, 1);
    lods_ratio_vec = zeros(num_architectures, 1);
    for i=1:num_architectures
        N_vec(i) = good_architectures(i).N;
        lods_ratio_vec(i) = max(good_architectures(i).L_i_full);
        architecture_str_vec{i} = good_architectures(i).arch;
    end
    pos_ratio_inds = find(lods_ratio_vec > 0);
    N_vec = N_vec(pos_ratio_inds);
    lods_ratio_vec = lods_ratio_vec(pos_ratio_inds);
    architecture_str_vec = architecture_str_vec(pos_ratio_inds);
    architecture_str_vec_short = architecture_str_vec;
    for i=1:length(architecture_str_vec_short)
        stop_ind = strfind(architecture_str_vec_short{i}, '(');
        if(~isempty(stop_ind))
            architecture_str_vec_short{i} = ...
                architecture_str_vec_short{i}(1:stop_ind-1);
        end
    end
    unique_arch = unique(architecture_str_vec_short);
    
    figure; hold on; % Plot minimal lods ratio as a function of N
    for i=1:length(unique_arch)
        arch_inds = ... % intersect(pos_ratio_inds, ...
            strmatch(unique_arch{i}, architecture_str_vec_short);
        cur_vec_to_plot = sortrows([N_vec(arch_inds) min(2, lods_ratio_vec(arch_inds))]); % sort rows
        unique_plot_min_inds = [1 find(diff(cur_vec_to_plot(:,1)))'+1];
        plot(cur_vec_to_plot(unique_plot_min_inds,1), ...
            cur_vec_to_plot(unique_plot_min_inds,2), color_vec(i), 'linewidth', 2); % dont plot anything too high
        plot(cur_vec_to_plot(unique_plot_min_inds,1), ...
            cur_vec_to_plot(unique_plot_min_inds,2), ['x' color_vec(i)], 'linewidth', 2); % dont plot anything too high
        
        %     plot(N_vec(arch_inds), min(2, lods_ratio_vec(arch_inds)), color_vec(i)); % dont plot anything to high
    end
    legend(mat_into_vec([unique_arch' ...
        empty_cell_to_empty_str(cell(length(unique_arch),1))]')', ...
        'location', 'southwest');
    title('minimal lods-ratio achieved for different architectures', ...
        'fontsize', font_size, 'fontweight', 'bold');
    xlabel('N loci', 'fontsize', font_size, 'fontweight', 'bold');
    ylabel('lods-ratio', 'fontsize', font_size, 'fontweight', 'bold');
    if(save_flag)
        my_saveas(gcf, [good_architectures_plot_file ...
            '_minimal_lods_ratio'], format_fig_vec);
    end
    
    figure; hold on; % Plot power for all architecture in a plot similar to Altchuler et al. for comparison
    plot_power_ind  = 1;
    for plot_power = [0.1 0.5 0.9]
        all_power_vec = []; % zeros(length(unique_arch), 1);
        all_maf_vec = []; % zeros(length(unique_arch), 1);
        all_sample_size_vec = []; % zeros(length(unique_arch), 1);
        all_lods_ratio_vec = [];
        params_struct.plot_power_vec = plot_power;
        subplot(2,2, plot_power_ind); hold on; plot_power_ind = plot_power_ind + 1;
        for i=1:length(good_architectures)
            good_power_inds = find( abs(good_architectures(i).power_marginal - ...
                params_struct.plot_power_vec) < 0.02 ); % find all values close in power
            if(~isempty(good_power_inds))
                all_power_vec = [all_power_vec ...
                    repmat( params_struct.plot_power_vec, 1, length(good_power_inds))];
                all_maf_vec = [all_maf_vec ...
                    repmat( good_architectures(i).MAF, 1, length(good_power_inds))];
                all_sample_size_vec = [all_sample_size_vec ...
                    params_struct.n_samples_vec(good_power_inds)];
                all_lods_ratio_vec = [all_lods_ratio_vec ...
                    repmat( max(good_architectures(i).L_i_full), 1, length(good_power_inds))];
            end
        end
        maf_ctr=1;
        for diff_mafs = unique(all_maf_vec)
            maf_inds = find(all_maf_vec == diff_mafs);
            plot(log(all_lods_ratio_vec(maf_inds)), log(all_sample_size_vec(maf_inds)), ...
                [color_vec(maf_ctr) '*']);
            maf_ctr = maf_ctr + 1;
        end
        x_ticks = num2str(exp(str2num(get(gca, 'xticklabels'))), 3)
        set(gca, 'xticklabels', x_ticks)
        y_ticks = num2str(round(exp(str2num(get(gca, 'yticklabels')))))
        set(gca, 'yticklabels', y_ticks)
        legend_vec = num2str_cell(num2cell(unique(vec2column(all_maf_vec))));
        for j=1:length(legend_vec)
            legend_vec{j} = ['MAF=' legend_vec{j}];
        end
        legend(legend_vec, 'location', 'southwest', 'fontsize', 8);
        xlabel('lods-ratio'); ylabel('sample size');
        title(['# samp. vs. lods-ratio for diff. MAF. Power = ' ...
            num2str(params_struct.plot_power_vec)]);
    end % end power plot (similar to Altchuler et al.)
    if(save_flag)
        my_saveas(gcf, [good_architectures_plot_file ...
            '_sample_size_power_vs_lods_ratio_power'], format_fig_vec); % ...
        %            num2str(params_struct.plot_power_vec)], format_fig_vec);
    end
    
    h_vec = zeros(num_architectures, 1); % Plot h on liability scale vs. disease scale
    h_add_vec = h_vec; h_liability_vec = h_vec;
    for i=1:num_architectures
        h_vec(i) = good_architectures(i).h_add;
        h_add_vec(i) = good_architectures(i).h_add;
         h_liability_vec(i) = good_architectures(i).h_liability;       
    end
    h_liability_transformed = h_add_vec; % CHANGE THIS FORMILAT !!! 
    figure; hold on; plot(h_add_vec, h_liability_vec, '*'); plot(0:0.1:1, 0:0.1:1, 'r');
    title('Heritability on liability vs. binary scales'); xlabel('Disease (binary) scale');
    ylabel('Liability scale');
    
    ctr=1;
    if(params_struct.plot_flag == 1) % only if flag says to plot
        for i=1:length(unique_arch)  % Make architecture-specific plots
            arch_inds = strmatch(unique_arch{i}, architecture_str_vec_short);
            
            for j=1:length(arch_inds) % loop on architectures with same ???
                N = length(good_architectures(arch_inds(j)).L_i_full);
                [cur_good_architectures_file cur_good_architectures_plot_file] = ...
                    get_architecture_file_names( ...
                    strdiff(strdiff(dir_from_file_name(good_architectures_plot_file), 'figures\'), ...
                    'figures/'), [arch_name_to_plot_name(unique_arch{i}) '_' params_struct.disease_type_str], ...
                    good_architectures(arch_inds(j)).N);
                
                figure;
                if(~seperate_plots)
                    subplot(2,2,1); hold on;   % Plot power for each architecture
                else
                    figure; hold on;
                end
                plot(params_struct.n_samples_vec, good_architectures(arch_inds(j)).power_marginal, ...
                    'linewidth', 3);
                if(isfield(good_architectures(arch_inds(j)), 'power_epistasis'))
                    plot(params_struct.n_samples_vec, ...
                        good_architectures(arch_inds(j)).power_epistasis, 'r--', ...
                        'linewidth', 2); % make -- line to still see the marginal
                    legend({['single-locus \alpha=' num2str(params_struct.power_alpha)], ...
                        ['pairwise-interaction \alpha=' ...
                        num2str(params_struct.power_pairwise_alpha)]}, ...
                        'location', 'northwest', 'fontsize', legend_font_size); %, ...
                else
                    legend('marginal-strongest-locus', 'location', 'northwest'); %, ...
                end
                title(['Det. power ' architecture_str_vec{arch_inds(j)} ' arch. (N=' num2str(N) ...
                    ') GWAS']); %, \alpha=' num2str(params_struct.power_alpha) ' (single), ' ...
                %                num2str(params_struct.power_pairwise_alpha) ' (pairwise)'], 'fontsize', font_size, 'fontweight', 'bold');
                xlabel('n (case+control)', 'fontsize', font_size, 'fontweight', 'bold');
                ylabel('power', 'fontsize', font_size, 'fontweight', 'bold');
                %             if(save_flag)
                %                 my_saveas(gcf, ...
                %                     [cur_good_architectures_plot_file '_power_arch_' num2str(arch_inds(j))], ...
                %                     format_fig_vec);
                %             end
                %figure;
                if(~seperate_plots)
                    subplot(2,2,2); hold on; % Plot prob. of getting disease for k one genotypes for each architecture
                else
                    figure; hold on;
                end
                errorbar(0:N, good_architectures(arch_inds(j)).mu_given_k_ones(N,:), ...
                    good_architectures(arch_inds(j)).mu_given_k_ones_std(N,:));
                plot(0:N, good_architectures(arch_inds(j)).mu_given_k_ones(N,:), 'r', 'linewidth', 2);
                plot(0:N, good_architectures(arch_inds(j)).mu_given_k_ones(N,:), 'r*');
                plot(0:N, repmat(good_architectures(arch_inds(j)).freq, N+1, 1), 'k--', 'linewidth', 2);
                legend('', 'Pr(dis.|k )', '', 'prevalence');
                xlabel('k loci', 'fontsize', font_size, 'fontweight', 'bold');
                ylabel('Prob. (disease)', 'fontsize', font_size, 'fontweight', 'bold');
                title(['Prob. disease, k of ' num2str(N) ' risk loci, ' ...
                    architecture_str_vec{arch_inds(j)} ' arch '], 'fontsize', font_size, 'fontweight', 'bold');
                axis([-0.1 N+0.1 -0.05 max(good_architectures(arch_inds(j)).mu_given_k_ones(N,:)) + 0.05]);
                %             if(save_flag)
                %                 my_saveas(gcf, [cur_good_architectures_plot_file ...
                %                     '_prob_disease_given_k_ones_arch_' num2str(arch_inds(j))], ...
                %                     format_fig_vec);
                %             end
                %figure;
                if(~seperate_plots)
                    subplot(2,2,3); hold on; % Plot simple binomial probability of risk alleles
                else
                    figure; hold on;
                end
                binom_probs = bernoulli_sum_prob(good_architectures(arch_inds(j)).f_vec);
                bar(0:N, binom_probs); title('Prob. of having k risk alleles');
                xlabel('k'); ylabel('Prob.');
                axis([-1 N 0 max(binom_probs) .* 1.1]);
                %             if(save_flag)
                %                 my_saveas(gcf, [cur_good_architectures_plot_file ...
                %                     '_prob_k_risk_alleles_arch_' num2str(arch_inds(j))], ...
                %                     format_fig_vec);
                %             end
                % figure;
                if(~seperate_plots)
                    subplot(2,2,4); hold on; % Plot relative risk
                else
                    figure; hold on;
                end
                max_generations = length(good_architectures(arch_inds(j)).relative_risk)
                plot(0:max_generations-1, good_architectures(arch_inds(j)).relative_risk, 'linewidth', 2);
                plot(0:max_generations-1, good_architectures(arch_inds(j)).relative_risk, '*');
                plot(0:max_generations-1, good_architectures(arch_inds(j)).freq*ones(max_generations,1), 'k--', ...
                    'linewidth', 2); % plot disease freq.
                plot(0:max_generations-1, good_architectures(arch_inds(j)).penet*ones(max_generations,1), 'r--', ...
                    'linewidth', 2); % plot disease freq.
                title('Risk to relative given disease for a person');
                ylabel('Prob.'); xlabel('relative degree');
                legend({'relative freq.', '', 'pop. freq.', 'penetrance'}, ...
                    'fontsize', legend_font_size);
                if(save_flag)
                    my_saveas(gcf, [cur_good_architectures_plot_file ...
                        '_all_four_figs_' num2str(good_architectures(arch_inds(j)).index)], ...
                        format_fig_vec);
                end
                
                figure; % new figure. Repeat the risk of k-of-N and binomial prob. when we know part of loci
                %                M_vec = [round(N/3) round(N/2)]; % try different # of known alleles
                M_vec = [round(N) round(N/2)]; % try different # of known alleles
                
                for m = 1:length(M_vec)
                    M = M_vec(m);
                    n_gwas = find(good_architectures(arch_inds(j)).power_marginal > M / N, 1);
                    if(isempty(n_gwas))
                        n_gwas = length(good_architectures(arch_inds(j)).power_marginal);
                    end
                    n_gwas = params_struct.n_samples_vec(n_gwas);
                    
                    
                    binom_probs = bernoulli_sum_prob(good_architectures(arch_inds(j)).f_vec(1:M));
                    if(m == 2) % here do the observed binomial - truncate edges
                       alpha_quantile = 0.01; 
                       cum_binom_probs = cumsum(binom_probs); 
                       x_binom = weighted_rand(binom_probs, 5000); 
                       max_M = max(x_binom); min_M = min(x_binom); % find boundaries 
                       binom_probs = hist(x_binom, min_M:max_M);
                    else
                        min_M = 0; max_M = M; % keep M as is                        
                    end
                    
                    subplot(2,2, m); hold on; % plot risk prob. for k ones out of M
                    %                     errorbar(0:M, good_architectures(arch_inds(j)).mu_given_k_ones(M,1:M+1), ...
                    %                         good_architectures(arch_inds(j)).mu_given_k_ones_std(M,1:M+1));
                    plot(min_M:max_M, good_architectures(arch_inds(j)).mu_given_k_ones(M,min_M+1:max_M+1), 'r', 'linewidth', 2);
                    plot(min_M:max_M, good_architectures(arch_inds(j)).mu_given_k_ones(M,min_M+1:max_M+1), 'r*');
                    plot(min_M:max_M, repmat(good_architectures(arch_inds(j)).freq, max_M-min_M+1, 1), 'k--');
                    legend('Pr(dis.|k ones)', '', 'prevalence');                    
                    xlabel('k loci', 'fontsize', font_size, 'fontweight', 'bold');
                    ylabel('Prob. (disease)', 'fontsize', font_size, 'fontweight', 'bold');
                    title(['Pr.(dis.), k of ' num2str(M) ' known loci, (GWAS n=' ...
                        num2str(n_gwas) ') ' architecture_str_vec{arch_inds(j)} ' arch '], 'fontsize', font_size, 'fontweight', 'bold');
                    axis([min_M-0.1 max_M+0.1 -0.05 max(good_architectures(arch_inds(j)).mu_given_k_ones(M,1:M+1))+0.05]);
                    subplot(2,2,m+2); hold on; % plot prob. of k ones out of M
                    bar(min_M:max_M, binom_probs); title(['Prob. k risk alleles of known ' num2str(M)]);
                    xlabel('k'); ylabel('Prob.');
                    axis([-1 max_M min_M max(binom_probs) .* 1.1]);
                end
                if(save_flag)
                    my_saveas(gcf, [cur_good_architectures_plot_file ...
                        '_risk_given_k_ones_of_M_known_loci_' ...
                        num2str(good_architectures(arch_inds(j)).index)], ...
                        format_fig_vec);
                end
                figure; hold on; % new figure: risk for a single pathway
                if(isfield(good_architectures(arch_inds(j)), 'mu_given_k_one_pathway'))
                    plot(0:k,  good_architectures(arch_inds(j)).mu_given_k_one_pathway, 'r');
                    plot(0:k,  good_architectures(arch_inds(j)).mu_given_k_one_pathway, 'r*');
                    title(['Risk of one pathway (' num2str(k) ' allels)']);
                    xlabel('number of risk loci'); ylabel('risk');
                end
                if(save_flag)
                    my_saveas(gcf, [cur_good_architectures_plot_file ...
                        '_risk_given_k_one_pathway_' ...
                        num2str(good_architectures(arch_inds(j)).index)], ...
                        format_fig_vec);
                end
                
                figure; hold on; % new figure:  family graph plot
                family_size = length(good_architectures(arch_inds(j)).family_tree);
                node_shapes = repmat([0 1 1 0], 1, ceil(family_size/4)); % get female/male shapes (box/circle)
                node_shapes = node_shapes(1:family_size);
                graph_draw(good_architectures(arch_inds(j)).family_tree, ...
                    'node_labels', num2str_cell(num2cell(good_architectures(arch_inds(j)).family_risk ./ ...
                    good_architectures(arch_inds(j)).freq), precision), ...
                    'node_shapes', node_shapes); % plot \lambda_R for each family member
                title('Relative risk-odds: person (grey) diseased given family members diseased');
                if(save_flag)
                    my_saveas(gcf, [cur_good_architectures_plot_file ...
                        '_family_risk_' num2str(good_architectures(arch_inds(j)).index)], ...
                        format_fig_vec);
                end
                
                ctr=ctr+1;
                if(mod(ctr, 10) == 0) % close figures
                    close all;
                end
            end % loop on architectures with same ???
        end % loop on unique_arch
    end % if plot_flag == 1
end % if num_architectures > 0

