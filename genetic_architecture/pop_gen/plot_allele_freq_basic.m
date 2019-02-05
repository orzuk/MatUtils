%
% Loop on curves and plot (bottom-up) - generic function
function plot_allele_freq_basic(x_vec_cell, p_vec_cell, plot_params) 

plot_params  = set_default_plot_params(plot_params);
my_symbol_vec = {'--', '-'}; % flip ordering (set integer powers as solid lines)
y_lim = [0 0];

if(plot_params.new_fig)
    figure;  
end
set(0, 'defaultFigureRenderer', 'painters'); % set render to show dots 

num_curves = length(x_vec_cell); 
for i = 1:num_curves
    
    %    tmp_color_ind = mod_max(6-floor(mod(i_s, 10)/2), 5);
    tmp_color_ind = mod_max(ceil(i/2), ceil(num_curves/2));
%    [~, i_s2] = min(abs(abs(D.s_grid)-abs(s_vec(i_s)))); % find closest s to current s with SFS pre-computed 
    if(iscell(x_vec_cell))
        plot_x_vec = x_vec_cell{i} ./ x_vec_cell{i}(end);
    else
        plot_x_vec = x_vec_cell ./ x_vec_cell(end);
    end
    if(iscell(p_vec_cell))
        plot_p_vec = p_vec_cell{i};
    else
        plot_p_vec = p_vec_cell(i,:);
    end
    if(plot_params.weighted) % weight by allele frequency
        if(i == 1)
            plot_params.ylabel_str = [plot_params.ylabel_str ' (weighted)'];
        end
        plot_p_vec = plot_p_vec .* plot_x_vec;
    end
    if(plot_params.normalize) % normalize distirbution
        plot_p_vec = plot_p_vec ./ sum(plot_p_vec);
    end
    if(plot_params.cum) % plot cumulative
        if(i == 1)
            plot_params.ylabel_str = [plot_params.ylabel_str ' (cum.)'];
        end
        if(plot_params.hist)
            plot_p_vec = cumsum(plot_p_vec); % cumsum_hist(plot_x_vec, plot_p_vec); % take histogram accounting for bins sizes. Need different normalization !!!
            if(plot_params.normalize)
                plot_p_vec = plot_p_vec ./ plot_p_vec(end);
            end
        else
            plot_p_vec = cumsum(plot_p_vec);
        end
        plot_params.legend_loc = 'southeast';
    else
        plot_params.legend_loc = 'northeast';
    end
    y_lim(1) = min(y_lim(1), min(plot_p_vec)); y_lim(2) = max(y_lim(2), max(plot_p_vec));
    
    %   max_diff_should_be_negative = max(diff(plot_p_vec))
    if(plot_params.log(1)) % log x
        if(plot_params.log(2)) % log y
            loglog(plot_x_vec, plot_p_vec, 'color', plot_params.color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i,2)}, 'linewidth', 2); hold on;
        else
            semilogx(plot_x_vec, plot_p_vec, 'color', plot_params.color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i,2)}, 'linewidth', 2); hold on;
        end
    else % no log on x
        if(plot_params.log(2)) % log y
            semilogy(plot_x_vec, plot_p_vec, 'color', plot_params.color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i,2)}, 'linewidth', 2); hold on;
        else
            plot(plot_x_vec, plot_p_vec, 'color', plot_params.color_vec{tmp_color_ind}, ...
                'linestyle', my_symbol_vec{mod_max(i,2)}, 'linewidth', 2); hold on;
        end
    end % log x
end % loop on i 

add_faint_grid(0.5, 0);
[h_leg, h_l] = legend(plot_params.legend, 'location', plot_params.legend_loc, 'fontsize', plot_params.font_size-4); legend('boxoff'); % just legend
xlabel('f (allele. freq.)'); ylabel(plot_params.ylabel_str); 

