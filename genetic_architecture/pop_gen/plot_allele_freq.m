% Plot allele frequency distributions for different models
% Input:
%
function plot_allele_freq(s_vec, Demographic_models, plot_params) % N, two_side_flag, log_flag, cum_flag, scale_mode, weight_flag)

AssignGeneralConstants; AssignRVASConstants;
eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
num_populations = length(Demographic_models);
num_s = length(s_vec);
[num_h, num_w] = num_to_height_width(num_populations);

plot_params  = internal_set_default_params(plot_params);

figure;
switch plot_params.figure_type % ALWAYS add legend !!!
    case 1 % CAF
        plot_params.ylabel_str = 'Combined allele frequency f_s';
        save_file = 'CAF_different_populations';
    case 2 % median of median allele frequency
        plot_params.ylabel_str = 'Median allele frequency';
        save_file = 'Median_IAF_different_populations';
    case 3 % mean of median allele frequency
        plot_params.ylabel_str = 'Mean of median allele frequency f_s';
        save_file = 'MeanMedian_IAF_different_populations';
    case 4 % again median ?? for non-zeros
        plot_params.ylabel_str = 'Median allele frequency f_s';
        save_file = 'Median_nonzero_IAF_different_populations';
end % switch figure type
for i_pop = 1:num_populations % loop on populations
    if(~isempty(Demographic_models{i_pop}))
        subplot(num_h, num_w, i_pop);
        internal_plot_allele_freq(s_vec, Demographic_models{i_pop}, plot_params);
    end
end
if(isfield(plot_params, 'figs_dir')) % save plot
    my_saveas(gcf, fullfile(plot_params.figs_dir, save_file), {'epsc', 'pdf', 'jpg'}); % NEW: add .jpg for Robert
end



% Internal plot
function internal_plot_allele_freq(s_vec, D, plot_params)

AssignGeneralConstants;
num_s = length(s_vec);
if(num_s <= 10)
    selection_color_vec = {'k', 'b', 'g', orange, 'r'}; % replace yellow with orange
else
    selection_color_vec = {'k', 'c', 'b', 'g', orange, 'r'}; % replace yellow with orange
end    
my_symbol_vec = {'--', '-'}; % flip ordering (set integer powers as solid lines)

s_legend_vec = s_vec_to_legend(s_vec);

for i_s = 1:num_s
    %    tmp_color_ind = mod_max(6-floor(mod(i_s, 10)/2), 5);
    tmp_color_ind = mod_max(ceil(i_s/2), ceil(num_s/2));
    
    [~, i_s2] = min(abs(abs(D.s_grid)-abs(s_vec(i_s))));
    if(iscell(D.SFS.x_vec))
        plot_x_vec = D.SFS.x_vec{i_s2} ./ D.SFS.x_vec{i_s2}(end);
    else
        plot_x_vec = D.SFS.x_vec ./ D.SFS.x_vec(end);
    end
    if(iscell(D.SFS.p_vec))
        plot_p_vec = D.SFS.p_vec{i_s2};
    else
        plot_p_vec = D.SFS.p_vec(i_s2,:);
    end
    
    if(plot_params.weighted) % weight by allele frequency
        plot_p_vec = plot_p_vec .* plot_x_vec;
    end    
    if(plot_params.normalize) % normalize distirbution
        plot_p_vec = plot_p_vec ./ sum(plot_p_vec);
    end
    if(plot_params.cum) % plot cumulative
        if(plot_params.hist)
            plot_p_vec = cumsum_hist(plot_x_vec, plot_p_vec); % take histogram accounting for bins sizes. Need different normalization !!! 
            if(plot_params.normalize)
                plot_p_vec = plot_p_vec ./ plot_p_vec(end); 
            end
        else            
            plot_p_vec = cumsum(plot_p_vec);
        end
    end
    %   max_diff_should_be_negative = max(diff(plot_p_vec))
    if(plot_params.log)
        loglog(plot_x_vec, plot_p_vec, 'color', selection_color_vec{tmp_color_ind}, ...
            'linestyle', my_symbol_vec{mod_max(i_s,2)}, 'linewidth', 2); hold on;
    else
        semilogx(plot_x_vec, plot_p_vec, 'color', selection_color_vec{tmp_color_ind}, ...
            'linestyle', my_symbol_vec{mod_max(i_s,2)}, 'linewidth', 2); hold on;
    end
end % loop on i_s 
ylabel(plot_params.ylabel_str, 'fontsize', 14); xlabel('f'); % tmp
title(D.name, 'fontsize', 14); % need to conver to nice name later
%add_faint_grid(0.5);
h_leg = legend(s_legend_vec, 'location', 'eastoutside'); % just legend
xlim(plot_params.xlim); ylim([min(plot_p_vec)*0.99, max(plot_p_vec)*1.01]);

%set(gca,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]); %ylim([10^(-6) 1]);
% PRoblem: legends appear twice - need to re-arrange order of grids to make
% them appear only once !!


% Internal function for setting defaults: density un-weighted
function plot_params  = internal_set_default_params(plot_params)

if(~isfield(plot_params, 'normalize')) % normalize distribution
    plot_params.normalize = 0;
end
if(~isfield(plot_params, 'cum')) % cumulative
    plot_params.cum = 0;
end
if(~isfield(plot_params, 'log')) % plot log-log
    plot_params.log = 0;
end
if(~isfield(plot_params, 'weighted')) % plot log-log
    plot_params.weighted = 0;
end
if(~isfield(plot_params, 'xlim')) % plot lim
    plot_params.xlim = [10^(-4) 1];
end
if(~isfield(plot_params, 'hist')) % hold distribution as histgoram
    plot_params.hist = 0;
end







% % % if(old_plot)
% % %
% % %     legend_vec = [repmat('s= -10^{', length(s_vec), 1) num2str(log10(abs(s_vec')),3) ...
% % %         repmat('}', length(s_vec), 1)];
% % %     legend_vec = cellstr(legend_vec);
% % %     legend_vec = strrep_cell(legend_vec, ' ', '');
% % %     legend_vec{1} = 's= 0'; % fix s=0
% % %     x_vec = (0:(2*N)) ./ (2*N);
% % %
% % %     figure;
% % %     for i=1:length(s_vec);
% % %         f_vec{i} = allele_freq_cumulative(x_vec, s_vec(i), N, two_side_flag, scale_mode, weight_flag);
% % %         semilogx(x_vec, f_vec{i}, [eric_color_vec(ceil(i/2)) my_symbol_vec{mod_max(i+1,2)}], 'linewidth', 2); hold on;
% % %     end
% % %     %legend(legend_vec, 'location', 'best'); legend('boxoff');
% % %     legend(legend_vec, 'position', [0.76 0.09 0.16 0.4]); legend('boxoff');
% % %
% % %     %    legend('frac. null', 'frac mutations captured');
% % %     xlabel('Derived allele frequency f'); ylabel('Proportion of alleles below f, \Psi_s(f)');
% % %     ylim([0 1.01*max(max_cell(f_vec))]);
% % %     if(log_flag)
% % %         xlim([10^(-4) 2]);
% % %     end
% % % end
