% Plot allele frequency distributions for different models
% Input:
% s_vec - vector of selection coefficients
% Demographic_models - structure with demographic information
% plot_params - parameters for plotting
%
function plot_allele_freq(s_vec, Demographic_models, plot_params) % N, two_side_flag, log_flag, cum_flag, scale_mode, weight_flag)

AssignGeneralConstants; AssignRVASConstants;
num_populations = length(Demographic_models);
num_s = length(s_vec);
[num_h, num_w] = num_to_height_width(num_populations);

plot_params  = internal_set_default_params(plot_params);
if(plot_params.new_fig)
    figure;  
end
set(0, 'defaultFigureRenderer', 'painters'); % set render to show dots 

switch plot_params.figure_type % ALWAYS add legend !!!
    case 1 % CAF
        plot_params.ylabel_str = 'Combined allele frequency $f_s$';
        save_file = 'CAF_different_populations';
    case 2 % median of median allele frequency
        plot_params.ylabel_str = 'Median allele frequency';
        save_file = 'Median_IAF_different_populations';
    case 3 % mean of median allele frequency
        plot_params.ylabel_str = 'Mean of median allele frequency $f_s$';
        save_file = 'MeanMedian_IAF_different_populations';
    case 4 % again median ?? for non-zeros
        plot_params.ylabel_str = 'Median allele frequency $f_s$';
        save_file = 'Median_nonzero_IAF_different_populations';
end % switch figure type
for i_pop = 1:num_populations % loop on populations
    if(~isempty(Demographic_models{i_pop}))
        plot_params.legend = {};
        subplot(num_h, num_w, i_pop); [i_w, i_h] = ind2sub([num_w, num_h], i_pop); 
        if(i_h == num_h)
            plot_params.legend = [plot_params.legend 'xlabel'];
        end
        if((i_w == 1) && (i_h == ceil(num_h/2))) % num_w)
            plot_params.legend = [plot_params.legend 'ylabel'];
        end
        if((i_w == num_w) && (i_h == ceil(num_h/2)))
            plot_params.legend = [plot_params.legend 'legend'];
        end
        internal_plot_allele_freq(s_vec, Demographic_models{i_pop}, plot_params);
    end
end
%orient landscape;

if(isfield(plot_params, 'figs_dir')) % save plot
    %orient landscape;    
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
y_lim = [0 0];

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
        if(i_s == 1)
            plot_params.ylabel_str = [plot_params.ylabel_str ' (weighted)'];
        end
        plot_p_vec = plot_p_vec .* plot_x_vec;
    end
    if(plot_params.normalize) % normalize distirbution
        plot_p_vec = plot_p_vec ./ sum(plot_p_vec);
    end
    if(plot_params.cum) % plot cumulative
        if(i_s == 1)
            plot_params.ylabel_str = [plot_params.ylabel_str ' (cum.)'];
        end
        if(plot_params.hist)
            plot_p_vec = cumsum_hist(plot_x_vec, plot_p_vec); % take histogram accounting for bins sizes. Need different normalization !!!
            if(plot_params.normalize)
                plot_p_vec = plot_p_vec ./ plot_p_vec(end);
            end
        else
            plot_p_vec = cumsum(plot_p_vec);
        end
    end
    y_lim(1) = min(y_lim(1), min(plot_p_vec)); y_lim(2) = max(y_lim(2), max(plot_p_vec));
    
    %   max_diff_should_be_negative = max(diff(plot_p_vec))
    if(plot_params.log(1)) % log x
        if(plot_params.log(2)) % log y
            loglog(plot_x_vec, plot_p_vec, 'color', selection_color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i_s,2)}, 'linewidth', 2); hold on;
        else
            semilogx(plot_x_vec, plot_p_vec, 'color', selection_color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i_s,2)}, 'linewidth', 2); hold on;
        end
    else % no log on x
        if(plot_params.log(2)) % log y
            semilogy(plot_x_vec, plot_p_vec, 'color', selection_color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i_s,2)}, 'linewidth', 2); hold on;
        else
            plot(plot_x_vec, plot_p_vec, 'color', selection_color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i_s,2)}, 'linewidth', 2); hold on;
        end
    end % log x
end % loop on i_s

xlim(plot_params.xlim); ylim([y_lim(1)*0.99, y_lim(2)*1.01]);
set(gca, 'XTick', logspace(log10(plot_params.xlim(1)), log10(plot_params.xlim(2)), log10(plot_params.xlim(2)) - log10(plot_params.xlim(1)) +1)); % change ticks
add_faint_grid(0.5, 0);


if(strmatch('ylabel', plot_params.legend))
    ylabel(plot_params.ylabel_str, 'fontsize', plot_params.font_size, 'interpreter', 'latex'); 
end
if(strmatch('xlabel', plot_params.legend))
    xlabel('$f$', 'interpreter', 'latex'); % tmp
end
title(strdiff(D.name, 'Fitted.'), 'fontsize', plot_params.font_size, 'interpreter', 'latex'); % need to conver to nice name later
if(strmatch('legend', plot_params.legend))
    [h_leg, h_l] = legend(s_legend_vec, 'location', 'eastoutside', 'fontsize', plot_params.font_size-4); legend('boxoff'); % just legend
    h_line=findobj(h_l,'type','line'); 
    lineXData = get(h_line, 'XData');  
    for j=1:2:length(lineXData) 
        lineXData{j}(2) = 0.65; lineXData{j}(1) = 0.3; 
        set(h_line(j), 'XData', lineXData{j});
    end
    % set(h_line, 'XData', lineXData);    
 %   set(h_leg, 'xcolor', [0.8 0.8 0.8], 'ycolor', [0.8 0.8 0.8]); % make legend invisible    
    pos_l = get(h_leg, 'position'); set(h_leg, 'position', [pos_l(1)+0.11 pos_l(2)-0.015 pos_l(3) pos_l(4)]);
end

% Internal function for setting defaults: density un-weighted
function plot_params  = internal_set_default_params(plot_params)
if(~isfield(plot_params, 'new_fig')) % normalize distribution
    plot_params.new_fig = 1;
end
if(~isfield(plot_params, 'normalize')) % normalize distribution
    plot_params.normalize = 0;
end
if(~isfield(plot_params, 'cum')) % cumulative
    plot_params.cum = 0;
end
if(~isfield(plot_params, 'log')) % default: plot semilogx
    plot_params.log = [1 0];
end
if(isscalar(plot_params.log))
    plot_params.log = [plot_params.log plot_params.log];
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
if(~isfield(plot_params, 'font_size')) % set default font size 
    plot_params.font_size = 14;
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
