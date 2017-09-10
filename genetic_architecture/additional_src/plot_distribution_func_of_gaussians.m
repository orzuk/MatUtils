% Plot density distribution of maximum of i.i.d. Gaussians
function plot_distribution_func_of_gaussians(N_vec, k_of_N_hist, k_of_N_bins_loc, fig_outfile)
AssignGeneralConstants;
num_N = length(N_vec);

for fig_ind = 1:2
    full_figure;
    
    for j=1:num_N % Just plot one of the heritability options
        if(fig_ind == 2)
            fig_str = '_normalized';
            [cur_p cur_x] = normalize_hist(k_of_N_bins_loc{1}{j}, k_of_N_hist{1}{j}, 1);
        else
            fig_str = '';
            cur_x = k_of_N_bins_loc{1}{j}; cur_p = k_of_N_hist{1}{j};
        end
        cur_p = smooth(cur_x, cur_p); % make figs look a bit nicer         
        plot(cur_x, cur_p, color_vec(j), 'linewidth', 2);
    end
    x_vec = -5:0.01:5; 
    legend_vec = vec2column(num2str_cell(num2cell(N_vec))); 
    if(fig_ind == 2)
        plot(x_vec, normpdf(x_vec), 'k--', 'linewidth', 4); % plot standard Gaussian 
        legend_vec{end+1} = 'Standard Gaussian'; 
    end
    xlabel('z'); ylabel('freq.'); legend(legend_vec);     
    xlim([-5 5]); 
    if(exist('fig_outfile', 'var') && (~isempty(fig_outfile)))
        my_saveas(gcf, [fig_outfile fig_str], format_fig_vec);
    end
end

