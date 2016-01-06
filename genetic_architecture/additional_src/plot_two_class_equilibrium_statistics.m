% Plot several statistics for the two-class model
%
% Input:
% two_class_stat_struct - structure with parameters
% w_x_null_mat - ??
% w_x_harmless - ??
% w_all - ??
% f_rare_vec - vector of allele frequencies f*
% frac_null_by_freq_cumulative - Matrix of Prob(allele null | freq. <= f). Notation:  \alpha_{s^*}(<= f)
% xi_prob_leq_f - Matrix of probability of allele being below prob. f. Notation: \xi( <= f)
% N - effective population size
% show_s_null - vector of specific values of s
% show_s_null_ind - vector of indices for specific values of s
% s_null_vec - vector of s for the null alleles
% rare_cumulative_per_gene - number of total rare variants per gene (1)
% alpha_vec - proportion of variants which are null
% title_str - string describing the two-class model parameters
% figs_dir - where to save figures
% new_figs_dir - where to save other figures
% figs_for_paper_flag - if this is 'on', just plot the figures needed for paper (if 'off', plot additional figures)
% demographic_models_struct - NEW! WE allow here more population models !!!
%
% Output:
% Compute several statistics for distribution of allele frequencies
%
function plot_two_class_equilibrium_statistics(two_class_stat_struct, w_x_null_mat, w_x_harmless, w_all, ...
    f_rare_vec, frac_null_by_freq_cumulative, xi_prob_leq_f, ...
    N, show_s_null, show_s_null_ind, s_null_vec, rare_cumulative_per_gene, alpha_vec, prevalence, ...
    title_str, figs_dir, new_figs_dir, plot_bayes_factor, figs_for_paper_flag, demographic_models_struct)

AssignGeneralConstants;
eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow) 
my_symbol_vec = {'-', '--'}; 


if(~exist('plot_bayes_factor', 'var') || isempty(plot_bayes_factor))
    plot_bayes_factor = 0; % plot independent of mixture coefficient
end
if(~exist('figs_for_paper_flag', 'var') || isempty(figs_for_paper_flag))
    figs_for_paper_flag = 0;
end
if(~exist('prevalence', 'var') || isempty(prevalence))
    prevalence=0.05; % assume a fixed prevalence (5%)
end
if(~exist('demographic_models_struct', 'var'))
    demographic_models_struct = []; demographic_models_struct.model_str = [];
end
demographic_models_struct.model_str{end+1} = 'equilibrium';
demographic_models_struct.num_models = length(demographic_models_struct.model_str);

%    legend_vec = [repmat('s^*= -', length(show_s_null), 1) num2str(show_s_null',3)];
legend_vec = [repmat('s= 10^{', length(show_s_null), 1) num2str(log10(show_s_null'),3) ...
    repmat('}', length(show_s_null), 1)];
legend_vec = cellstr(legend_vec);
legend_vec = strrep_cell(legend_vec, ' ', '');
legend_vec{1} = 's= 0'; % fix s=0

prob_null_given_x_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,end);


if(~figs_for_paper_flag) % plot many figures
    for log_x_flag = 0:1 % Plot frac null and frac captured on the same plot
        figure;
        tmp_title_str = [strrep(title_str, ...
            'Detect. power rare', 'Frac. null') ' ; frac. null (solid); frac. captured (dashed)'];
        tmp_title_str = strdiff(tmp_title_str, 'cum. freq.');
        %    title(tmp_title_str);
        switch log_x_flag
            case 0 % linear plot
                x_log_str = '';
                plot(f_rare_vec, frac_null_by_freq_cumulative{1}(show_s_null_ind,:), 'linewidth', 2); hold on;
                plot(f_rare_vec, xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene, '--', 'linewidth', 2);
            case 1 % logarithmic plot
                x_log_str = '_log';
                semilogx(f_rare_vec, frac_null_by_freq_cumulative{1}(show_s_null_ind,:), 'linewidth', 2); hold on;
                semilogx(f_rare_vec, xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene, '--', 'linewidth', 2);
        end
        
        
        
        legend(legend_vec, 4-log_x_flag);   legend('boxoff');
        %    legend('frac. null', 'frac mutations captured');
        xlabel('Max allele freq. f^*'); ylabel('frac. null (solid) ; frac. captured (dashed)');
        ylim([0 1.01*max(max(frac_null_by_freq_cumulative{1}(show_s_null_ind,:)))]);
        my_saveas(gcf, fullfile(figs_dir, 'power', ...
            ['rare_variants_null_fraction_as_function_of_freq_cutoff_alpha_' num2str(alpha_vec,2) x_log_str]), {'epsc', 'pdf'});
    end % loop on log flag
end % if ~figs_for_paper

close all;

for i_d = 1:demographic_models_struct.num_models    % New!!! plot also for other populations !!!!
    switch demographic_models_struct.model_str{i_d}
        case 'equilibrium'
            plot_f_vec = f_rare_vec;
        otherwise
            demographic_models_struct.data{i_d} = load(demographic_models_struct.file_names{i_d}); % this gives both s-vec and cumulatives
            plot_f_vec = demographic_models_struct.data{i_d}.x_vec;
    end
    for plot_type = 0:1 % 0 - frac. which are null. 1 - frac. which are captured
        switch plot_type
            case 0  % fraction which are null
                if(plot_bayes_factor)
                    y_vec = w_x_null_mat{1}(show_s_null_ind,:) ./ ...
                        repmat(w_x_harmless{1}, length(show_s_null_ind), 1); %   y_vec ./ alpha_vec;
                    y_str = 'relative depletion for nulls';
                else
                    switch demographic_models_struct.model_str{i_d}
                        case 'equilibrium'
                            y_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:);
                        otherwise
                            y_vec = demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) ./ ...
                                ( demographic_models_struct.data{i_d}.p_vec(end:-1:1,:) + ...
                                repmat(demographic_models_struct.data{i_d}.p_vec(end,:), 10, 1) ); % This should go DOWN!!!
                    end
                    y_str = 'frac. null';
                end
                cur_file_str = 'frac_null_given_freq_less_than_cutoff_alpha_'
                cur_title_str = ['Frac. of null alleles for alleles with freq. below f^*. \alpha=' num2str(alpha_vec)];
                %                prob_null_given_x_vec = y_vec(:,end);
            case 1  %-cumulative fraction captured
                switch demographic_models_struct.model_str{i_d}
                    case 'equilibrium'
                        y_vec = xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene;
                    otherwise
                        y_vec = bsxfun(@plus, demographic_models_struct.data{i_d}.p_vec(end:-1:1,:), ...
                            demographic_models_struct.data{i_d}.p_vec(end,:));
                end % switch
                cur_file_str = 'rare_variants_freq_less_than_cutoff_alpha_'
                cur_title_str = ['Frac. of all alleles with freq. below f^*. \alpha=' num2str(alpha_vec)];
                y_str = 'frac. below';
        end % switch plot type
        cur_title_str = [cur_title_str ' . ' demographic_models_struct.model_str{i_d}];
        for log_x_flag = 1 % 0:1
            for log_y_flag = 0 % 0:1
                %%            if( (figs_for_paper_flag && log_x_flag) ) %  && (~log_y_flag) )
                figure;     % New plot separately and on log scale
                %%            else
                %%                subplot(2,2, log_x_flag + 2*log_y_flag+1);
                %%            end
                switch log_x_flag
                    case 0
                        x_log_str = '';
                        switch log_y_flag
                            case 0
                                y_log_str = '';
                                %    plot(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2); hold on;
                                plot(f_rare_vec, y_vec, 'linewidth', 2);
                            case 1
                                y_log_str = ' (log)';
                                %    semilogy(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2);  hold on;
                                semilogy(f_rare_vec,y_vec, 'linewidth', 2);
                        end
                    case 1
                        x_log_str = ' (log)';
                        switch log_y_flag
                            case 0
                                y_log_str = '';
                                %    semilogx(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2);  hold on;
                                for j=1:size(y_vec, 1)
                                    semilogx(plot_f_vec, y_vec(j,:), [eric_color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}],  'linewidth', 3); hold on;
                                end
                            case 1
                                y_log_str = ' (log)';
                                %    loglog(f_rare_vec, frac_null_by_freq_cumulative(show_s_null_ind,:), 'linewidth', 2);  hold on;
                                for j=1:size(y_vec, 1)
                                    loglog(plot_f_vec, y_vec, [eric_color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}], 'linewidth', 3); hold on; 
                                end
                        end
                end % switch log x
                xlabel(['Max allele freq. f^*' ]); % x_log_str
                ylabel([y_str y_log_str]);
                if(log_x_flag)
                    xlim([10^(-4) 1]);
                end
                if(log_y_flag) % set the same axis for x and y
                    ylim([10^(-4) 1]);
                end
                %    ylim([0 1.01*max(max(frac_null_by_freq_cumulative(show_s_null_ind,:)))]);
                if(figs_for_paper_flag)
                    switch plot_type
                        case 0
                            output_fig_file = 'figure1b';
                            title('Fig 1b', 'fontsize', 14, 'fontweight', 'b');
                            if(plot_bayes_factor)
                                ylim([0 1.01]);
                            end
                        case 1 % plot combined null_neutral freq.
                            if(~log_y_flag)
                                ylim([0 1]);
                            end
                            output_fig_file = ['figure1c' y_log_str];
                            title('Fig 1c', 'fontsize', 14, 'fontweight', 'b');
                    end
                else
                    switch log_y_flag
                        case 0
                            if(plot_type == 0)
                                output_fig_file = fullfile(demographic_models_struct.model_str{i_d}, ...
                                    ['fig2c_' cur_file_str num2str(alpha_vec,2)]);
                            else
                                output_fig_file = fullfile(demographic_models_struct.model_str{i_d}, ...
                                    ['fig2b_' cur_file_str num2str(alpha_vec,2)]);
                            end
                            
                        case 1
                            output_fig_file = fullfile('power', [cur_file_str num2str(alpha_vec,2)]);
                    end % switch plot_type
                    suptitle(cur_title_str);
                end
                %            if(~figs_for_paper_flag)
                legend(legend_vec,3);  legend('boxoff'); % +plot_type*(~figs_for_paper_flag));
                %            end
                if(plot_type == 1)
                    ylim([0 1]);
                else
                    ylim([0 alpha_vec .* 1.05]); 
                end
                add_faint_grid();
                
                my_saveas(gcf, fullfile(new_figs_dir, output_fig_file), {'epsc', 'pdf'});
            end % loop on log-y
        end % loop on log-x
    end % loop on plot type
end % loop on population models

% end % optional plots

close all;


for i_d = 1:demographic_models_struct.num_models
    for log_x_flag = 1 % 0:1 % Plot null and harmless alleles separately
        for neutral_flag = 0 % 0:1 % 0: null; 1: neutral
            switch neutral_flag
                case 0
                    neutral_str = 'null';
                    cur_frac_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:);
                case 1
                    neutral_str = 'neutral';
                    cur_frac_vec = 1-frac_null_by_freq_cumulative{1}(show_s_null_ind,:);
            end
            switch demographic_models_struct.model_str{i_d}
                case 'equilibrium'
                    plot_f_vec = f_rare_vec;
                    cur_y_vec = cur_frac_vec .* ...
                        xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene;
                otherwise
                    plot_f_vec = demographic_models_struct.data{i_d}.x_vec;
                    cur_y_vec = demographic_models_struct.data{i_d}.p_vec(end:-1:1,:);
            end
            for normalize_flag = 1 % 0:1 % choose if to normalize (divide by
                switch normalize_flag
                    case 0 % do nothing
                        normalize_str = '';
                    case 1
                        normalize_str = 'normalized_';
                        cur_y_vec = cur_y_vec ./ repmat(cur_y_vec(:,end), 1, size(cur_y_vec, 2)); % normalized to get to one
                end
                figure;    % one plot
                tmp_title_str = [strrep(title_str, ...
                    'Detect. power rare', 'Frac. null captured') ' ;'];
                tmp_title_str = strdiff(tmp_title_str, 'cum. freq.');
                %    title(tmp_title_str);
                switch log_x_flag
                    case 0 % linear plot
                        x_log_str = '';
                        plot(plot_f_vec, cur_y_vec, 'linewidth', 2);
                        %             plot(f_rare_vec, frac_null_by_freq_cumulative{1}(show_s_null_ind,:), 'linewidth', 2); hold on;
                        %             plot(f_rare_vec, xi_prob_leq_f{1}(show_s_null_ind,:) ./ rare_cumulative_per_gene, '--', 'linewidth', 2);
                    case 1 % logarithmic plot
                        x_log_str = '_log';
                        for j=1:size(cur_y_vec, 1)
                            semilogx(plot_f_vec, cur_y_vec(j,:), [eric_color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}], 'linewidth', 3); hold on;
                        end
                end
                %            legend_vec = [repmat('s^*= -', length(show_s_null), 1) num2str(show_s_null',3)];
                legend(legend_vec, 3-log_x_flag); legend('boxoff');
                %    legend('frac. null', 'frac mutations captured');
                xlabel('Max allele freq. f^*'); ylabel('frac. captured');
                ylim([0 1.01*max(max(cur_y_vec))]);
                if(log_x_flag)
                    xlim([10^(-4) 1]);
                end
                if(figs_for_paper_flag)
                    if(neutral_flag)
                        title('Fig 1b', 'fontsize', 14, 'fontweight', 'b');
                    else
                        ylabel('frac. alleles');
                        title('Fig 1a', 'fontsize', 14, 'fontweight', 'b');
                    end
                    output_fig_file = 'figure1a';
                else
                    %                 output_fig_file = fullfile('power', ['rare_variants_' normalize_str neutral_str '_fraction_captured_as_function_of_freq_cutoff_alpha_' ...
                    %                     num2str(alpha_vec,2) x_log_str]); % fig doesn't depend on alpha
                    output_fig_file = fullfile('new_eric', demographic_models_struct.model_str{i_d}, ...
                        ['fig2a_rare_variants_' normalize_str neutral_str '_cumulative_distribution_f_x' ...
                        x_log_str]);
                end
                add_faint_grid();
                title(tmp_title_str); 
                my_saveas(gcf, fullfile(figs_dir, output_fig_file), {'epsc', 'pdf'});
            end % loop on normalization flag
        end % loop on neutral flag
    end % loop on log flag
end % loop on models

close all; 

% New: Plot cumulative variance explained
%figure; cum_het_vec =

% semilogx(f_rare_vec, cur_y_vec, 'linewidth', 2);
if(figs_for_paper_flag)
    var_expl_vec = zeros(length(show_s_null), length(f_rare_vec));
    for i=1:length(show_s_null)
        run_i = i
        for j=1:length(f_rare_vec)
            %    var_expl_vec(i,:) = XXX(f_rare_vec, g =
            %   var_expl_vec(i,:) =  exp(allele_freq_spectrum(f_rare_vec, -show_s_null(i), N, 0, 'log', 1));
            var_expl_vec(i,j)=  absorption_time_by_selection(show_s_null(i), 1, N, 0, f_rare_vec(j), -2);
            %        absorption_time_by_selection(s, theta, N, x_min, x_max, weight_flag)
        end
        var_expl_vec(i,:) = var_expl_vec(i,:) ./ var_expl_vec(i,end); % normalize
    end
    % for i=1:length(show_s_null)
    %     var_expl_vec(i,:) =
    figure;
    %plot(f_rare_vec, var_expl_vec, 'linewidth', 2);
    semilogx(f_rare_vec, var_expl_vec, 'linewidth', 2);
    xlim([10^(-4) 1]);
    xlabel('Max allele freq. f^*'); ylabel('cum. var. expl.');
    % legend(legend_vec, 3);
    title('Fig 1d', 'fontsize', 14, 'fontweight', 'b');
    my_saveas(gcf, fullfile(figs_dir, 'figure1d'), {'epsc', 'pdf'});
end % if plot figs for paper



if(~figs_for_paper_flag)
    % New! Plot variance explained by one locus as a function of s
    %tmp_num_loci_vec = -S_vec .* integral_phi_x' ./ (-S_vec - (1-exp(S_vec)))
    tmp_var_explained_per_gene_vec = ... % This is the var explained per GENE
        allele_freq_cumulative(1-1/(2*N), -s_null_vec, N, 0, 'linear', 'var') - ...
        allele_freq_cumulative(1/(2*N), -s_null_vec, N, 0, 'linear', 'var'); % Start with variance explained per gene
    
    tmp_var_explained_per_polymorphic_locus_vec = ... % This is the var explained per GENE
        ( allele_freq_cumulative(1-1/(2*N), -s_null_vec, N, 0, 'linear', 'var') - ...
        allele_freq_cumulative(1/(2*N), -s_null_vec, N, 0, 'linear', 'var') ) ./ ...
        ( allele_freq_cumulative(1-1/(2*N), -s_null_vec, N, 0, 'linear', 0) - ...
        allele_freq_cumulative(1/(2*N), -s_null_vec, N, 0, 'linear', 0) ); % Start with variance explained per gene
    tmp_num_genes_vec = 1 ./ tmp_var_explained_per_gene_vec;
    tmp_num_polymorphic_loci_vec = 1 ./ tmp_var_explained_per_polymorphic_locus_vec;
    
    
    figure; loglog(s_null_vec, tmp_var_explained_per_gene_vec, 'linewidth', 2 ); hold on;
    loglog(s_null_vec, tmp_var_explained_per_polymorphic_locus_vec, 'r', 'linewidth', 2 );
    xlabel('-s'); ylabel('Var. Explained');
    legend('gene', 'polymorphic allele');
    my_saveas(gcf, fullfile(figs_dir, 'power', ...
        'var_explained_per_gene_as_function_of_selection'), ...
        {'epsc', 'pdf'});
    figure; loglog(s_null_vec, tmp_num_genes_vec, 'linewidth', 2 ); hold on;
    loglog(s_null_vec, tmp_num_polymorphic_loci_vec, 'r', 'linewidth', 2 );
    xlabel('-s'); ylabel('Num. Genes');
    legend({'gene', 'polymorphic allele'},4);
    my_saveas(gcf, fullfile(figs_dir, 'power', ...
        'num_genes_needed_as_function_of_selection'), ...
        {'epsc', 'pdf'});
    
    % figure; loglog(s_null_vec, tmp_var_explained_per_polymorphic_locus_vec, 'linewidth', 2 );
    % xlabel('-s'); ylabel('Var. Explained');
    % my_saveas(gcf, fullfile(figs_dir, 'power', ...
    %     'rare_variants_var_explained_per_polymorphic_allele_as_function_of_selection'), ...
    %     {'epsc', 'pdf'});
    % figure; loglog(s_null_vec, tmp_num_polymorphic_loci_vec, 'linewidth', 2 );
    % xlabel('-s'); ylabel('Num. Genes');
    % my_saveas(gcf, fullfile(figs_dir, 'power', ...
    %     'num_polymorphic_alleles_needed_as_function_of_selection'), ...
    %     {'epsc', 'pdf'});
    
    
    
    % Alternative: Numeric version
    x_vec = (1:N-1) ./ (N);
    s = -0.0001; y_vec = allele_freq_spectrum(x_vec, s, N, 0);
    z_vec = normalize_hist(x_vec, y_vec);
    VV = integral_hist(x_vec, x_vec .* (1-x_vec) .* z_vec), 1 / VV
    
    two_class_stat_struct.median_x = max(two_class_stat_struct.median_x, 0); % can't have negative values
    for plot_flag = 0:1 % what to plot
        for log_flag = 0:1
            figure;
            switch log_flag
                case 0
                    log_str = '';
                    switch plot_flag
                        case 0
                            plot(s_null_vec, two_class_stat_struct.mean_x, 'linewidth', 2); hold on;
                            plot(s_null_vec, two_class_stat_struct.median_x, 'r', 'linewidth', 2);
                        case 1
                            plot(s_null_vec, two_class_stat_struct.mean_het_x, 'g', 'linewidth', 2);
                    end
                case 1
                    log_str = '_log';
                    switch plot_flag
                        case 0
                            semilogx(s_null_vec, two_class_stat_struct.mean_x, 'linewidth', 2); hold on;
                            semilogx(s_null_vec, two_class_stat_struct.median_x, 'r', 'linewidth', 2);
                        case 1
                            semilogx(s_null_vec, two_class_stat_struct.mean_het_x, 'g', 'linewidth', 2);
                    end
            end
            %title('Mean allele freq. as function of s');
            switch plot_flag
                case 0
                    het_str = 'allele_freq';
                    legend('mean-allele-freq.', 'median-allele-freq.');
                case 1
                    het_str = 'heterozygosity';
                    legend('mean-heterozygosity');
            end
            xlabel('-s'); % ylabel('mean-allele-freq.');
            my_saveas(gcf, fullfile(new_figs_dir, 'power', ...
                ['mean_' het_str '_given_s_' num2str(alpha_vec,2) log_str] ), {'epsc', 'pdf'});
        end % loop on log-flag
    end % loop on plot flag
    
    R = [ show_s_null'  vec2column(two_class_stat_struct.mean_x(show_s_null_ind)) ...
        vec2column(two_class_stat_struct.median_x(show_s_null_ind)) vec2column(two_class_stat_struct.mean_het_x(show_s_null_ind)) prob_null_given_x_vec]; % Save also a table
    R = [{'S', 'mean-allele-freq', 'median-allele-freq.', 'mean-heterozygosity', ['prob-null (\alpha=' num2str(alpha_vec) ')']}' num2str_cell(num2cell(R'),2)]';
    savecellfile(R, fullfile(new_figs_dir, 'power', 'mean_allele_freq_tab.txt'));
    
    
    
    % Plot mean number of alleles present in a single chromosome
    mu=1.6*10^(-8); L=625; %% 1500; % take typical parameters for mutation rate and gene length. We want to get mu*L=10^(-5)
    
    
    figure; models_legend_vec = []; model_ctr=1; % New! plot in this figure all populations together
    for i_d = 1:demographic_models_struct.num_models    % New!!! plot also for other populations !!!!
        switch demographic_models_struct.model_str{i_d}
            case 'equilibrium'
                num_alleles_per_chr_vec = 2 .* N .* mu .* L .* 2 .* ...
                    two_class_stat_struct.normalization_factor_x, 'b', 'linewidth', 2;
                plot_s_null_vec = s_null_vec;
            case {'varsel1', 'varsel2', 'varselection1', 'varselection2'}
                continue; % we don't plot for variable selection
            otherwise
                plot_s_null_vec = max(10^(-6), demographic_models_struct.data{i_d}.s_vec); % add a point at 10^(-6) to represent zero
                num_alleles_per_chr_vec = demographic_models_struct.data{i_d}.num_alleles_per_chrom;
        end
        semilogx(plot_s_null_vec, num_alleles_per_chr_vec,  color_vec(i_d), 'linewidth', 2); hold on;
        models_legend_vec{model_ctr} = demographic_models_struct.model_str{i_d}; model_ctr=model_ctr+1;
    end % loop on demographic model
    xlabel('s'); ylabel('# alleles-per-chrom.');
    title(['Expected num. derived alleles per chrom. $\mu_{G,null}=' num2str(mu*L) '$. ' ...
        demographic_models_struct.model_str{i_d}], 'interpreter', 'latex');
    add_faint_grid(); legend(models_legend_vec); legend('boxoff'); 
    my_saveas(gcf, fullfile(new_figs_dir, 'all_models', ...
        'fig1a_mean_num_alleles_per_chromosome_as_function_of_s'), {'epsc', 'pdf'});
    
    
    
    
    %end % if test power
    
    % Plot power (what is proxy for power?)
    y_vec = frac_null_by_freq_cumulative{1}(show_s_null_ind,:).^2 .* xi_prob_leq_f{1}(show_s_null_ind,:);
    y_vec = y_vec ./ repmat(max(y_vec,[],2), 1, size(y_vec, 2))
    figure; plot(f_rare_vec, y_vec, 'linewidth', 2);
    legend(legend_vec,4); title('Power proxy = (frac. null.)^2 * frac-kept, as function of allele-freq. cutoff');
    xlabel('Max allele freq. f^*'); ylabel('Proxy for power');
    my_saveas(gcf, fullfile(figs_dir, 'power', ...
        ['rare_variants_power_vs_freq_cutoff_alpha_'  num2str(alpha_vec,2)]), {'epsc', 'pdf'});
    
    
    % Plot isolines for variance explained with effect size and s on the two axis
    GRR_vec = 1:0.01:8; % genetic relative risk
    f_vec = 0.0000001; % temp (should change f, but the resulting beta isn't very sensitive to allele freq.
    total_desired_var_explained = logspace(-4,-1, 7);
    
    
    
    
    for i_d = 1:demographic_models_struct.num_models    % New!!! plot also for other populations !!!!
        switch demographic_models_struct.model_str{i_d}
            case 'equilibrium'
                mean_het_vec_per_gene = 2.*N.*mu .* L .* 2 .* two_class_stat_struct.mean_het_x .* two_class_stat_struct.normalization_factor_x;
                %                num_alleles_per_chr_vec = 2 .* N .* mu .* L .* 2 .* two_class_stat_struct.normalization_factor_x, 'b', 'linewidth', 2;
                plot_s_null_vec = s_null_vec;
            otherwise
                plot_s_null_vec = demographic_models_struct.data{i_d}.s_vec; % only few val
                demographic_models_struct.data{i_d}.density_p_vec = zeros(size(demographic_models_struct.data{i_d}.p_vec));
                demographic_models_struct.data{i_d}.mean_het_x = zeros(length(plot_s_null_vec), 1);
                for i_s = 1:length(plot_s_null_vec)
                    demographic_models_struct.data{i_d}.density_p_vec(i_s,:) = [vec2row(diff(demographic_models_struct.data{i_d}.p_vec(i_s,:))) 0];
                    demographic_models_struct.data{i_d}.mean_het_x(i_s) = ...
                        demographic_models_struct.data{i_d}.num_alleles_per_chrom(i_s) * ...
                        mean_hist( demographic_models_struct.data{i_d}.x_vec', demographic_models_struct.data{i_d}.density_p_vec(i_s,:));
                end
                mean_het_vec_per_gene = 2 .* (demographic_models_struct.data{i_d}.num_alleles_per_chrom - ...
                    demographic_models_struct.data{i_d}.mean_het_x); % demographic_models_struct.data{i_d}.num_alleles_per_chrom;
        end
        GRR_isoline_vec = zeros(length(plot_s_null_vec), length(total_desired_var_explained));
        for j=1:length(total_desired_var_explained)
            beta_vec = sqrt(total_desired_var_explained(j) ./ mean_het_vec_per_gene); % var. explained proportional to beta^2
            GRR_isoline_vec(:,j) = beta_to_genetic_relative_risk(beta_vec, f_vec, prevalence);
        end
        isocurve_legend_vec = cellstr( [repmat('$V=10^{', length(total_desired_var_explained), 1) num2str(log10(total_desired_var_explained')) ...
            repmat('}$', length(total_desired_var_explained), 1)])
        
        
        figure;
        for j=1:size(GRR_isoline_vec,2)
            semilogx(plot_s_null_vec, GRR_isoline_vec(:,j), [color_vec(ceil(j/2)) my_symbol_vec{mod_max(j,2)}],  'linewidth', 3); hold on;
        end
        ylim([0.999 10]);
        legend(isocurve_legend_vec, 'interpreter', 'latex', 'location',  'NorthWest');  legend('boxoff');
        xlabel('s'); ylabel('GRR'); title(['Isocurves for variance explained per gene. $\pi=' num2str(prevalence) '$. ' ...
            demographic_models_struct.model_str{i_d}], 'interpreter', 'latex');
        add_faint_grid();
        my_saveas(gcf, fullfile(new_figs_dir, demographic_models_struct.model_str{i_d}, ...
            'fig1b_var_explained_per_gene_as_function_of_s_and_GRR_isocurves'), {'epsc', 'pdf'});
    end % loop on demographic models
end % if figs_for_paper

%%%%%%%%%%%%%%%%%%%%%% Figure 2: Plot things as function of s, the selection coefficient %%%%%%%%%%%%%%%%%%%%%%





close all;



%%%%%%%%%%%
% Make some simple plots of density
if(~figs_for_paper_flag)
    cumulative_psi_s = cell(3,1);
    f_rare_vec(end) = 0.999999999999;
    for plot_type=2 % loop on different
        ctr=1;
        for j=show_s_null_ind
            cumulative_psi_s{plot_type}(ctr,:) = phi_s_integral(f_rare_vec,-s_null_vec(j)*4*N, plot_type-1);
            cumulative_psi_s{plot_type}(ctr,:) = ...
                cumulative_psi_s{plot_type}(ctr,:) - cumulative_psi_s{plot_type}(ctr,1);
            cumulative_psi_s{plot_type}(ctr,:) = ...
                cumulative_psi_s{plot_type}(ctr,:) ./ cumulative_psi_s{plot_type}(ctr,end); % normalize to one
            ctr=ctr+1;
            switch plot_type
                case 1 % here do \int_{0}^x psi_s(x)
                    plot_str = '\phi_s(x)';
                    
                case 2 % here do \int_0^x x*psi_s(x)
                    plot_str = 'x\phi_s(x)';
                    
                case 3 % here do \int_0^x x^2*psi_s(x)
                    plot_str = 'x^2\phi_s(x)';
                    
            end % switch
        end % loop on selection coefficient s
        file_str = strrep(plot_str, ' ', '_');
        file_str = strrep(file_str, '\', '_');
        file_str = strrep(file_str, '^', '_');
        file_str = strrep(file_str, '(', '_');
        file_str = strrep(file_str, ')', '_');
        
        figure; % Plot ..
        for log_x = 0:1
            for log_y = 0:1
                log_y_str = '';
                subplot(2,2,log_x+log_y*2+1);
                switch log_x
                    case 0
                        log_x_str = '';
                        switch log_y
                            case 0
                                plot(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                                
                            case 1
                                log_y_str = '(log)';
                                semilogy(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                        end
                    case 1
                        log_x_str = '(log)';
                        switch log_y
                            case 0
                                log_x_str = '(log)';
                                semilogx(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                                
                            case 1
                                log_y_str = '(log)';
                                loglog(f_rare_vec, cumulative_psi_s{plot_type}, 'linewidth', 2);
                                legend(legend_vec, 4);
                        end
                end
                title(['Cumulative fraction of distribution $\int_{0}^{f^*} ' plot_str ' dx$'], 'interpreter', 'latex');
                xlabel(['Allele Freq. f^* ' log_x_str]); ylabel(['Cumulative density ' log_y_str]);
            end % log x
        end % log y
        my_saveas(gcf, fullfile(new_figs_dir, 'power', ...
            ['cumulative_allele_freq_' file_str]), {'epsc', 'pdf'});
    end % loop on plot type
end % optional plots

close all;
