% Plot power for MLT model
%
% Input:
% isoheritability_flag - whether we force heritability to be the same for all models
% N - # of pathways/liabilities
% star_h_x_vals - values we indicate in stars (should match the main text)
%
function n_samples_power_table = ...
    plot_heritability_power_figures( isoheritability_flag, N, star_h_x_vals)
AssignGeneralConstants;
format_fig_vec = 'epsc';

power_outfile = '../../common_disease_model/figs/epistasis_power/MLT_detection_power_different_h_x';
%fig_outdir = '../../common_disease_model/docs/figs/new'; % New! figs for paper
fig_outdir = '../../common_disease_model/docs/pnas/genetic_interactions/figs'; % New! figs for paper

%load(['power_MLT_' num2str(N) '_1_data_iso_' num2str(isoheritability_flag) ]); % pathway ?
load_paper_fig_data = 0;
if(load_paper_fig_data)
    X=load('pathway_epistasis_gof_power_MLT_3_1_data_iso_0.mat'); % pathway??
    Y=load('Big_power_MLT_3_1_data_iso_0_marginal_and_pairwise_epistasis.mat'); % pathway
    load('Big_power_MLT_3_1_data_iso_0_marginal_and_pairwise_epistasis.mat'); % pathway
    try_test_inds = union(X.try_test_inds, Y.try_test_inds);
    power_pathway_epistasis_cell = Y.power_pathway_epistasis_cell;
    power_pathway_epistasis_cell{10} = X.power_pathway_epistasis_cell{10};
else
    %    load(['power_MLT_data_iso_1_' num2str(isoheritability_flag)]);
    load(['power_LP_disease_' num2str(N) '_1_data_iso_' num2str(isoheritability_flag) ]);
end


num_n_samples = length(n_samples_vec);
power_table = power_pathway_epistasis_cell; % cell(num_stats,1);
power_stat_str = str2title(test_name_vec);
alpha_vec = power_pathway_alpha_vec; % [power_pathway_alpha power_pathway_alpha power_pathway_alpha power_pathway_alpha power_pathway_alpha];

power_vec = [0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99]; % desired powers for iso-power curves
%power_half_ind = [1 3 5 7 9]; % (power_vec == 0.5); % index for 50% power
power_half_ind = [3 5 7]; % (power_vec == 0.5); % index for 50% power


%power_vec = [0.05 0.1 0.25 0.5 0.75 0.9 0.95]; % desired powers for iso-power curves

num_powers = length(power_vec);
num_h = length(h_x_one_locus_vec);
good_inds = [1 2 5:size(power_table,1)]; % New: good inds exclude some of the statistics
good_inds = intersect(good_inds, try_test_inds); % Take only tests for which statistics were computed
power_table = power_table(good_inds,:); % take only a few test statistics. Next few are
NCP_table = non_centrality_parameter(good_inds,:)
power_stat_str = power_stat_str(good_inds);
alpha_vec = alpha_vec(good_inds);
num_stats = size(power_table,1); % length(power_stat_str);
num_mu = length(mu_vec); % number of different prevalences

%1:num_h; % [1:2:9 10]; % pick only some of the h_x's
n_samples_power_table = zeros(num_powers, num_stats,num_h,num_mu);
R = [];
for mu_ind = 1:num_mu % loop on different prevalences
    mu = mu_vec(mu_ind);
    for k=1:num_h % loop on different heritabilities
        for i=1:num_powers % Fill table. loop on power cutoffs
            for j=1:num_stats % loop on different power statistics
                n_ind = find(power_table{j,mu_ind}(k,:) > power_vec(i), 1);
                if(~isempty(n_ind)) % find the sample size giving a pre-specified power
                    if(n_ind == 1) % here even the smallest n has enough power
                        n_samples_power_table(i,j,k,mu_ind) = n_samples_vec(n_ind);
                    else % linear inerpolation between two values of n
                        alpha = (power_vec(i) - power_table{j,mu_ind}(k,n_ind)) / ...
                            (power_table{j,mu_ind}(k,n_ind-1) - power_table{j,mu_ind}(k,n_ind));
                        n_samples_power_table(i,j,k,mu_ind) = round( ...  % Take linear interpolation
                            alpha*n_samples_vec(n_ind-1) +  (1-alpha)*n_samples_vec(n_ind) );
                    end
                else
                    n_samples_power_table(i,j,k,mu_ind) = max(n_samples_vec); % Take the maximum value (we knot the true value is larger)
                end
                cur_NCP_vec(j) = NCP_table{j,mu_ind}(k,end) / n_samples_vec(end);
            end % loop on statistic used
        end % loop on power
        [max_inds_I max_inds_J] = find(n_samples_power_table(:,:,k,mu_ind) == max(n_samples_vec)); % here we just know that n>something
        
        n_samples_power_table_str = num2str_cell(num2cell(n_samples_power_table(:,:,k,mu_ind)));
        for i=1:length(max_inds_I)
            n_samples_power_table_str{max_inds_I(i),max_inds_J(i)} = ...
                ['>' n_samples_power_table_str{max_inds_I(i),max_inds_J(i)}];
        end
        n_samples_power_table_str = [num2str_cell(num2cell(alpha_vec))' ...
            num2str_cell(num2cell(cur_NCP_vec))' n_samples_power_table_str'];
        n_samples_power_table_str = [power_stat_str' n_samples_power_table_str];
        GRR = heritability_to_genetic_relative_risk(... % diploid 
            h_x_one_locus_vec(k)/2, 'liability', freq, mu);
        n_samples_power_table_str = [[['h_x=' num2str(100*h_x_one_locus_vec(k),2) '%, (RAF=' ...
            num2str(freq*100,3) '%, GRR=' num2str(GRR,3) ')'] ...
            cell(1,num_powers+2)]' n_samples_power_table_str']';
        if(k==1) % first heritability
            n_samples_power_table_str = [[{'Test:', '\alpha:', 'NCP'} cell(1,num_powers)]' n_samples_power_table_str']';
            n_samples_power_table_str = [[{'', ... % ['h_x=' num2str(100*h_x_one_locus_vec(k),2) '%']
                '', 'Power:'} num2str_cell(num2cell(100*power_vec),[],[],1)]' ...
                n_samples_power_table_str']';
        end
        n_samples_power_table_str = [cell(1,num_powers+3)' n_samples_power_table_str']';
        if(k==1)
            n_samples_power_table_str = [[['Sample size needed for detection in balanced case-control study. ' ...
                'Model: LP(N=' num2str(N) ',K=1,h=100%,\mu=' num2str(100*mu) '%). RAF=' num2str(100*freq,3) '%'] cell(1,num_powers+2)]' n_samples_power_table_str']';
            if(mu_ind > 1)
                n_samples_power_table_str = [cell(2,num_powers+2)' n_samples_power_table_str']';
            end
        end
        R = [R' n_samples_power_table_str']';
    end % loop on heritability
    %    R = [R' ' '
end % end loop on prevalence

savecellfile(R, [power_outfile '.txt'], [], 1); % Save power table
h_z=h_x_one_locus_vec; % this is the heritability on liability scale assuming ONE liability!!!


for mu_ind = 1:num_mu % New: Plot iso-power curves: what sample sizes give 50% power for each statistic
    full_figure(0);
    legend_vec = cell(num_stats,1);
    mu = mu_vec(mu_ind);
    if(num_stats == 3)
        j_vec = [2 1 3];
    else
        j_vec = 1:num_stats;
    end
    for j=j_vec % [2 1 3] % :num_stats % :-1:1 % loop on different power stats (we should have here both 'pathway', single-locus and pairwise loci
        if(isempty(strfind_cell(test_name_vec(good_inds), 'pathway')))
            cur_symbol = symbol_vec{(num_stats-j)*2+1}; % (2-j)*2+1
        else
            cur_symbol = symbol_vec{mod_max(j-1,3)}; %   symbol_vec{j}; % symbol_vec{2*j-1};
        end
        
        %        for quant_ind = power_half_ind
        legend_vec = num2str_cell(num2cell(100*power_vec(power_half_ind)),[],[],1);
        %        legend_vec{j} = [power_stat_str{j} ' (\alpha=' num2str(alpha_vec(j), 3) ')'];
        y_vec = reshape(n_samples_power_table(power_half_ind,j,:,mu_ind), ...
            length(power_half_ind), num_h)' + 1;
        semilogy(100*h_x_one_locus_vec,  y_vec, ...
            cur_symbol, 'linewidth', 2); % color_vec(j),
        hold on;
        if(j <= length(star_h_x_vals)) % plot stars corresponding to examples in the text
            switch j
                case 1
                    cur_symbol = 'd';
                case 2
                    cur_symbol = '*';
                case 3
                    cur_symbol = 'o';
            end
            [match_h_val match_h_ind] = min(abs(100*h_x_one_locus_vec - star_h_x_vals(j))); % find where is the star
            plot(100*h_x_one_locus_vec(match_h_ind), ... % text example
                n_samples_power_table(5,j,match_h_ind,mu_ind), ...
                ['r' cur_symbol],  'linewidth', 2, 'markersize', 12); % plot star belonging to the text.
        end
        if(j==num_stats)
            hold on;
        end
        %        end
    end % loop on statistic used
    legend(legend_vec, 3); legend boxoff;
    %     title(['sample size reaching 50% power for different statistics. MLT(' ...
    %         num2str(N) ', 1, \mu=' num2str(mu*100,3) '%, h=100%)']);
    if(isempty(strfind_cell(test_name_vec(good_inds), 'pathway')))
        title_str = 'fig3a';
    else
        title_str = 'supp-fig5'; %  'fig3'; % this fig was moved to supp. info. 
    end
    include_stat_str = 0;
    if(include_stat_str)
        title([str2title(cell2vec(test_name_vec(good_inds), ', ')) '   ' ...
            str2title(cell2vec(symbol_vec(1:3), ', '))]);
    else
        title([repmat(' ', 1, 100) str2title(title_str)], ...
            'fontsize', 16, 'fontweight', 'bold');
    end
    xlabel('h_i^2 (%)');
    ylabel('n');   % sample size (cases+controls)
    xlim([0  max(h_x_one_locus_vec)*100]);
    y_ticks = get(gca, 'YTIck'); y_ticks(end)=y_ticks(end)*0.9989999;
    y_tick_labels = get(gca, 'YTIcklabel');
    y_tick_labels = [repmat('$10^', size(y_tick_labels,1), 1) ...
        y_tick_labels repmat('$', size(y_tick_labels,1), 1)]
    %     set(gca, 'YTicklabel', y_tick_labels)
    %     set(gca, 'YTick', y_ticks)
    ylim([min(n_samples_power_table(:)) max(n_samples_power_table(:))*1]);
    my_saveas(gcf, fullfile(fig_outdir, title_str), format_fig_vec);
    %     my_saveas(gcf, [power_outfile '_iso_power_curves_MLT_' num2str(N) ...
    %         '_1_mu_' num2str(mu,3)], format_fig_vec);
    
end % loop on mu


%return; % Just do the main plots for now .. .

for mu_ind = 1:num_mu % loop on different prevalences
    for k=num_h:num_h % 1:num_h % New: make a figure for each heritability
        cur_power_table = zeros(num_n_samples, num_stats);
        for i=1:num_n_samples
            for j=1:num_stats
                cur_power_table(i,j) = power_table{j,mu_ind}(k,i);
            end
        end
        full_figure(0);
        semilogx(n_samples_vec, 100*cur_power_table(:,1:min(7,length(good_inds))), ...
            'linewidth', 2); hold on;
        if(length(good_inds) > 7)
            hold on;
            semilogx(n_samples_vec, 100*cur_power_table(:,8:end), ...
                'linewidth', 2, 'linestyle', ':'); hold on;
        end
        xlabel('Sample size'); ylabel('detection power (%)');
        
        %    for i=1:length(h_z) % Compute effect sizes (main and interaction)
        compute_h_in_plot = k
        [z_mu_LT z_var_LT] = z_stats_given_two_x_MLT(1, 1, mu, h_z(k), ...
            h_x_one_locus_vec(k)); % replaced h_x_MLT_one_locus_vec(mu_ind,k)
        [z_mu_MLT z_var_MLT] = z_stats_given_two_x_MLT(N, K, mu, h_z(k), ...
            h_x_MLT_one_locus_vec(mu_ind,k));
        var_explained_main(k) = (z_mu_LT*(1-z_mu_LT) - z_var_LT) / z_var_LT;
        var_explained_interaction(k) = (z_var_LT - z_var_MLT) / z_var_LT;
        %    end
        
        
        %     legend_vec = num2str_cell(num2cell(h_x_one_locus_vec(good_inds)*100),2,[],1);
        %     for i=1:length(legend_vec)
        cur_title_vec = [', V_{main}=' num2str(100*var_explained_main(k),3) ...
            '%, V_{interaction}=' num2str(100*var_explained_interaction(k),3) '%'];
        %     end
        title(['Epistasis det. power MLT(' ...
            num2str(N) ',' num2str(K) ',\mu=' num2str(mu_vec(mu_ind)*100),'%,h=' num2str(h_x*100,3) '%) with' ...
            ' h_x^2=' num2str(100*h_x_one_locus_vec(k),3) '% (1 liability)' cur_title_vec]); % :-main, solid-epistasis. ' ...
        %        'x-epistasis no-interaction.']); %  v-epistasis model-based. ^-epistasis model-based-no-interaction.']);
        legend_vec = cell(num_stats,1);
        for j=1:num_stats
            legend_vec{j} = [power_stat_str{j} ' (\alpha=' num2str(alpha_vec(j), 3) ')'];
        end
        legend(legend_vec  , 2); % Show values for different heritabilities
        my_saveas(gcf, [power_outfile '_' num2str(h_x_one_locus_vec(k),3)], format_fig_vec);
    end % loop on heritabilities
end % loop on prevalence

return;

for j=1:num_stats  % loop on stats. Make one figure per statistic
    full_figure(0); % plot assuming # different h_x <= 7
    mu_title_vec = [];
    for mu_ind = 1:num_mu % plot all mu's together
        cur_power_table = zeros(num_n_samples, num_h);
        for i=1:num_n_samples
            for k=1:num_h
                cur_power_table(i,k) = power_table{j,mu_ind}(k,i);
            end
        end
        semilogx(n_samples_vec, 100*cur_power_table, symbol_vec{mu_ind}, 'linewidth', 2); hold on;
        %     if(length(good_inds) > 7)
        %         hold on;
        %         semilogx(n_samples_vec, 100*cur_power_table(:,8:end), ...
        %             'linewidth', 2, 'linestyle', ':'); hold on;
        %     end
        mu_title_vec = [mu_title_vec '\mu=' num2str(100*mu_vec(mu_ind),3) '% (' ...
            symbol_str_vec{mu_ind} ') '];
    end % loop on prevalence
    xlabel('Sample size'); ylabel('detection power (%)');
    title(['Epistasis det. power MLT(' ...
        num2str(N) ',' num2str(K) ',\mu' ',h=' num2str(h_x*100,3) '%). Stat: ' ...
        power_stat_str{j} ' (\alpha=' num2str(alpha_vec(j), 3) '). ' mu_title_vec]);
    legend_vec = cell(num_h,1);
    for k=1:num_h
        legend_vec{k} = ['h_x=' num2str(h_x_one_locus_vec(k)*100,3) '%'];
    end
    legend(legend_vec  , 2); % Show values for different heritabilities
    my_saveas(gcf, [power_outfile '_MLT_' num2str(N) '_1_mu_' num2str(mu) '_h_' num2str(h_x) ...
        '_stat_' power_stat_str{j} '_alpha_' num2str(alpha_vec(j), 3)], format_fig_vec);
    
end % loop on different statistics








% %     full_figure(0); % Plot power as function of heritability of pathway (many genes)
% %     semilogx(n_samples_vec, 100*[power_marginal power_epistasis ...
% %         power_pathway_epistasis_mat_gaussian_genotypes' ...
% %         power_pathway_epistasis_mat_gaussian_genotypes_goodness_of_fit'], 'linewidth', 2);
% %     %    plot(log10(n_samples_vec), power_epistasis, 'r');
% %     legend({'main-effect, \alpha=5*10^{-8}', ...
% %         ['epistasis, \alpha=' num2str(power_pairwise_alpha,2)], ...
% %         ['pathway-epistasis, \alpha=' num2str(power_pathway_alpha,2)], ...
% %         ['pathway-gaussian-epistasis, \alpha=' num2str(power_pathway_alpha,2)], ...
% %         ['pathway-gaussian-epistasis-goodness-of-fit, \alpha=' num2str(power_pathway_alpha,2)]}, 2);
% %     xlabel('Sample size (log_{10})'); ylabel('detection power (%)');
% %     title(['Detection power for the MLT(' ...
% %         num2str(N) ',' num2str(K) ',' num2str(mu*100),'%,' num2str(h_x*100,3) '%) model with ' ...
% %         num2str(num_loci) ' loci in each liability']);
% %     my_saveas(gcf, '../../common_disease_model/figs/MLT_detection_power', format_fig_vec);

if(exist('pathway_grr_vec', 'var'))
    full_figure; plot(n_samples_vec, pathway_grr_vec-1, '*');
    plot(n_samples_vec, h_x_one_locus*num_loci*power_marginal, 'r*');
    title('Total pathway GRR');
    legend('GRR-1', 'h_x^2'); xlabel('n_{samples}');
end
