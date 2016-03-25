% Plot var. explained distribution (quantiles) for different values of
% selection coefficient
% 
% Input:  
% quant_vec - vector of quantile values 
% N - effective population size
% figs_dir - where to save output figures 
%
function plot_var_explained_by_allele_freq_quantiles(quant_vec, N, figs_dir)

AssignGeneralConstants;
x_grid = 0.5.*(1:N)./N;
for s_grid_ind = 1:2 %
    if(s_grid_ind == 1)
        s_grid = -logspace(-6, 0, 100); % make a finer grid 
    else
        s_grid = -[0:0.00001:0.005];
    end
    num_s = length(s_grid);
    
    mean_freq = zeros(num_s,1);
    std_freq = zeros(num_s,1);
    median_freq = zeros(num_s,1);
    mean_var_explained_freq= zeros(num_s,1);
    std_var_explained_freq= zeros(num_s,1);
    median_var_explained_freq= zeros(num_s,1);
    num_quant = length(quant_vec);
    quantile_freq = zeros(num_s, num_quant);
    quantile_var_explained = zeros(num_s, num_quant);
    for i=1:length(s_grid) % can take a long time if we've got many s values
        y_grid = exp(allele_freq_spectrum(x_grid, s_grid(i), N, 1, 'log'));
        y_var_grid = y_grid .* x_grid .* (1-x_grid);
        mean_freq(i) = mean_hist(x_grid, y_grid);
        median_freq(i) = median_hist(x_grid, y_grid);
        std_freq(i) = std_hist(x_grid, y_grid);
        mean_var_explained_freq(i) = mean_hist(x_grid, y_var_grid);
        std_var_explained_freq(i) = std_hist(x_grid, y_var_grid);
        median_var_explained_freq(i) = median_hist(x_grid, y_var_grid);
        
        for j=1:length(quant_vec)
            quantile_freq(i,j) = quantile_hist(x_grid, y_grid, quant_vec(j)); % Get also quantiles:
            quantile_var_explained(i,j) = ...
                quantile_hist(x_grid, y_var_grid, quant_vec(j)); % Get also quantiles:
        end
        if(s_grid_ind == 1)  % make an error term
            quantile_freq(i,:) = abs(quantile_freq(i,:) - median_freq(i));
            quantile_var_explained(i,:) = ...
                abs(quantile_var_explained(i,:) - median_var_explained_freq(i));
        end
    end
    s_label = num2str_cell(num2cell(-s_grid), 2)'
    big_S_grid = 4*N.*s_grid;
    big_S_label = num2str_cell(num2cell(-big_S_grid), 2)'
    error_bar_plot_mode = 0;
    format_fig_vec = {'fig', 'jpg', 'epsc', 'pdf'};
    for plot_ind = s_grid_ind:3 % start with allele freq, then var. expl and then log
        full_figure(0); % hold on;
        ylabel('s');
        if(s_grid_ind == 1) % plot whiskers
            if(plot_ind == 1) % plot ...
                line_x = 0.01;
                xlabel('allele. freq. dist.');
                title('Allele Freq. Dist. for different s');
                plot_str = 'allele_freq_dist_by_s';
                %            herrorbar(mean_freq, 1:num_s,  std_freq);
                switch error_bar_plot_mode
                    case 1 % Just plot bars 
                        herrorbar(median_freq, 1:num_s,  ...
                            quantile_freq(:,1),  quantile_freq(:,end), 'r');
                        herrorbar(median_freq, 1:num_s,  ...
                            quantile_freq(:,2),  quantile_freq(:,end-1), 'g');
                        herrorbar(median_freq, 1:num_s,  ...
                            quantile_freq(:,3),  quantile_freq(:,end-2), 'm');
                        herrorbar(median_freq, 1:num_s,  ...
                            quantile_freq(:,4),  quantile_freq(:,end-3));
                    case 0 % here connect the lines
                        for j=1:4 % why do we go up to four here??
                            plot([median_freq quantile_freq(:,j) ...
                                quantile_freq(:,end-j+1)], 1:num_s, color_vec(j+4))
                        end
                end % switch how to plot (error bars or not)
                xlabel('MAF'); 
            else
                if(plot_ind == 2) % plot ... 
                    xlabel('MAF');
                    plot_x = median_var_explained_freq;
                    line_x = 0.01;
                    plot_str = 'var_explained_dist_by_s';
                else % at 3 plot on log scale
                    xlabel('MAF (log-scale)');
                    plot_x = log10(median_var_explained_freq);
                    line_x = log10(0.01);
                    plot_str = 'var_explained_dist_by_s_log';
                end
                title('Var. Expl. Dist. for different s');
                
                switch error_bar_plot_mode
                    case 1
                        %            herrorbar(mean_var_explained_freq, 1:num_s,  std_var_explained_freq);
                        herrorbar(plot_x, 1:num_s,  ...
                            quantile_var_explained(:,1),  quantile_var_explained(:,end), 'r');
                        herrorbar(plot_x, 1:num_s,  ...
                            quantile_var_explained(:,2),  quantile_var_explained(:,end-1), 'g');
                        herrorbar(plot_x, 1:num_s,  ...
                            quantile_var_explained(:,3),  quantile_var_explained(:,end-2), 'm');
                        herrorbar(plot_x, 1:num_s,  ...
                            quantile_var_explained(:,4),  quantile_var_explained(:,end-3));
                    case 0 % here connect the lines
                        for j=1:4
                            semilogy(plot_x, -s_grid, color_vec(j+4), 'linewidth', 2); % just to get the legends right
                            if(j == 1)
                                hold on;
                            end
                        end
                        semilogy(plot_x, -s_grid, 'k', 'linewidth', 4);
                        for j=1:4
                            semilogy([plot_x-quantile_var_explained(:,j) ...
                                plot_x+quantile_var_explained(:,end-j+1)], -s_grid, ...
                                color_vec(j+4), 'linewidth', 2);
                        end
                end % switch how to plot (error bars or not)
                xlabel('MAF');
            end
            line([line_x line_x], [10^(-6) 10^2], 'color', 'k', 'linestyle', '--'); % Mark the 'rare' boundary line
            
            switch error_bar_plot_mode % Still need to fix the legend here !!!!
                case 1
                    legend({'1%', '99%', '5%', '95%', '10%', '90%', '25%', '75%'});
                case 0
                    legend({'1%, 99%', '5%, 95%', '10%, 90%', '25%, 75%' 'median (50%)'});
            end
        else % here plot lines
            if(plot_ind == 2)
                xlabel('MAF');
                plot_x = median_var_explained_freq;
                line_x = 0.01;
                plot_str = 'var_explained_vs_s_plot';
                plot(quantile_var_explained, -s_grid); % one plot !!!
                %                 for j=1:num_quant
                %                     plot(quantile_var_explained(:,j), -s_grid, [color_vec(abs(j-5)+1)]);
                %                 end
            else
                xlabel('MAF (log-scale)');
                plot_x = log10(median_var_explained_freq);
                line_x = log10(0.01);
                plot_str = 'var_explained_vs_s_plot_log';
                semilogx(quantile_var_explained, -s_grid);
                %                 for j=1:num_quant
                %                     plot(log10(quantile_var_explained(:,j)), -s_grid, [color_vec(abs(j-5)+1)]);
                %                 end
            end
            legend({'1%', '5%', '10%', '25%', '50%', '75%', '90%', '95%', '99%', });
            line([line_x line_x], [0 log10(max(abs(s_grid)))], 'color', 'k', 'linestyle', '--'); % Mark the 'rare' boundary line
        end % if s_grid_ind
        
        if(plot_ind ~= 2) % for 2 just use semilog which fixes everything
            max_y_tick = get(gca, 'ytick');
            if(s_grid_ind == 1)
                res = 4; % show onlt 10^x where x is integer
                set(gca, 'ytick', 1:res:length(s_grid));
                s_scale = get(gca, 'ytick');
                set(gca, 'yticklabel', s_label);
                if(plot_ind == 3) % log-scale on x-axis
                    x_tick = get(gca, 'xtick');
                    x_tick_labels = get(gca, 'xticklabel');
                    new_x_tick_labels = num2str_cell(num2cell(10.^x_tick), 2)
                    %            new_x_tick_labels = [repmat('10^(', 10, 1) x_tick_labels repmat(')', 10, 1)];
                    set(gca, 'xticklabel', new_x_tick_labels);
                end
                add_right_yticks([max_y_tick(1) max_y_tick(end)], 'S=4Ns');
                set(gca, 'ytick', 1:length(s_grid));
                set(gca, 'yticklabel', big_S_label);
            else
                add_right_yticks([max_y_tick(1) max_y_tick(end)], 'S=4Ns');
                y_tick_labels = get(gca, 'yticklabel');
                big_S_label = num2str_cell(num2cell(abs(-4*N.*max_y_tick)), 3)'
                %            set(gca, 'ytick', 1:length(s_grid));
                set(gca, 'yticklabel', big_S_label);
            end
        else % here 2nd plot - copy y-axis to the right side
            y_lim = [10^(-6) 10^(-1)];
            ylim(y_lim);
            
            y_ticks = get(gca, 'YTick'); % Get right axis ticks % [-10^(-2) -10^(-6)]; % get(gca, 'YTick');
            y_tick_labels = get(gca, 'YTickLabel');
            log_scale = 1;
            if(log_scale)
                y_tick_labels = [repmat('10^{', length(y_tick_labels), 1) ...
                    y_tick_labels repmat('}', length(y_tick_labels), 1)];
            end
            axesPosition = get(gca,'Position');          %# Get the current axes position
            hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YLim',y_lim,...     %#   ... and a different scale
                'YAxisLocation','right',...  %#   ... located on the right
                'XTick',[],...               %#   ... with no x tick marks
                'YTick', [], ...
                'YTickLabel', [], ...
                'Box','off');                %#   ... and no surrounding box
            %    y_ticks = get(gca, 'YTick'); % Get right axis ticks % [-10^(-2) -10^(-6)]; % get(gca, 'YTick');
            if(log_scale)
                y_ticks = linspace(y_lim(1), y_lim(2), length(y_ticks));
            end
            my_yticklabels(hNewAxes, y_ticks,  mat2cell(y_tick_labels, ones(length(y_ticks),1)), ...
                'YAxisLocation', 'right');

        end
        my_saveas(gcf, fullfile(figs_dir, plot_str), format_fig_vec);
    end % loop on plot ind
end % loop on what type of s/slope and plot to perform


