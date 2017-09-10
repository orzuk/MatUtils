% Plot density distribution of maximum of i.i.d. Gaussians
% 
% Input: 
% N_vec - vector of N, number of gaussians 
% max_flag - plot distirbution of maximum (default) or minimum
% fig_outfile - save figures in this file
% 
function plot_max_of_gaussians(N_vec, max_flag, fig_outfile)

AssignGeneralConstants;
format_fig_vec = 'epsc';
add_kurtosis = 0; % how to display legend 

num_N = length(N_vec);
if(~exist('max_flag', 'var'))
    max_flag = 'MAX';
end

switch upper(max_flag)
    case 'MAX'
    case 'MIN'
end

x_vec = -100:0.001:100; % Generate a dense x vector
mu = zeros(num_N,1); sigma = zeros(num_N,1); k = zeros(num_N,1); s = zeros(num_N,1); % moments

for i=1:num_N % loop on number of gaussians
    y_vec{i} = N_vec(i) .* normcdf(x_vec).^(N_vec(i)-1) .* normpdf(x_vec);  % Generate distribution
    [y_vec_normalized{i} bin_locs{i}] = normalize_hist(x_vec, y_vec{i}, 1);
    k(i) = kurtosis_hist(bin_locs{i}, y_vec_normalized{i})% compute kurtosis
    k_unnormalized(i) = kurtosis_hist(x_vec, y_vec{i});
    s(i) = skewness_hist(bin_locs{i}, y_vec_normalized{i})% compute skewness
    s_unnormalized(i) = skewness_hist(x_vec, y_vec{i});
    mu(i) = mean_hist(x_vec, y_vec{i}); sigma(i) = sqrt(var_hist(x_vec, y_vec{i}));
end
gumbel_hist = evpdf(-x_vec); [gumbel_hist_normalized gumbel_bins_loc] = normalize_hist(x_vec, gumbel_hist, 1); % Plot limit distribution
k_unnormalized(end+1) = kurtosis_hist(x_vec, gumbel_hist);
k(end+1) = kurtosis_hist(gumbel_bins_loc, gumbel_hist_normalized); % compute kurtosis for Gumbel. Should be 2.4
s_unnormalized(end+1) = skewness_hist(x_vec, gumbel_hist);
s(end+1) = skewness_hist(gumbel_bins_loc, gumbel_hist_normalized); % compute skewness for Gumbel. Should be ??
k(abs(k) < 0.000001) = 0; s(abs(s) < 0.000001) = 0; % round small values
mu(1) = 0; mu(end+1) = 10^9; sigma(end+1) = 0; % the limiting (un-normalzied) distribution
for fig_ind = 1:4 % first: main, second: supp, third: un-normalized supp, fourth: scaled and un-scaled together
    full_figure; % Plot normalized trait distributions
    if(add_kurtosis)
        plot(0, 0, 'color', 'w'); % dummy plot for empty legend 
    end
    switch fig_ind
        case 1 % plot normalized, Main figure, remove inifnity
            x_mat = cell2vec(bin_locs(1:end), [], 0)';
            y_mat = cell2vec(y_vec_normalized(1:end), [], 0)';
            fig_str = '_main'; legend_place = 1;
        case 2 % plot normalized.
            x_mat = [cell2vec(bin_locs, [], 0)']; %  gumbel_bins_loc'];
            y_mat = [cell2vec(y_vec_normalized, [], 0)']; %  gumbel_hist_normalized'];
            fig_str = ''; legend_place = 1;
        case 3 % plot unnormalized
            x_mat = repmat(x_vec, length(y_vec), 1)';
            y_mat = cell2vec(y_vec, [], 0)';
            fig_str = '_unnormalized'; legend_place = 2;
    end
    switch upper(max_flag)
        case 'MAX'
        case 'MIN'
            x_mat = -x_mat;
    end

    
    %     for i=1:num_N % loop on number of gaussians
    %     if(i < 4)
    %         j=i;
    %     else
    %         j=i+1;
    %     end
    %    plot(bin_locs{i}, y_vec_normalized{i}, color_vec(mod_max(j, 96)), 'linewidth', 2); % plot distributions on top of each others. Avoid black
    %    plot(x_mat, y_mat, 'linewidth', 2); % plot distributions on top of each others. Avoid black
    plot(x_mat(:,1), y_mat(:,1), 'k', 'linewidth', 2); % plot standard gaussian (thicker?)
    plot(x_mat(:,2:end), y_mat(:,2:end), 'linewidth', 2); % plot again null distribution
    if(fig_ind == 2)
        plot(gumbel_bins_loc, gumbel_hist_normalized, 'k--', 'linewidth', 2); % plot again infinity Gumbel distribution
        y_mat = [y_mat gumbel_hist_normalized'];
    end
    legend_vec = [num2str_cell(num2cell(N_vec)) '\infty'];
    max_len = max(length_cell(legend_vec(1:size(x_mat,2))));
    
    if(add_kurtosis)
        for j=1:num_N+1 % add skewness and kurtosis to legend
            if(length(legend_vec{j}) < max_len)
                legend_vec{j} = [repmat(' ', 1, 2*(max_len-length(legend_vec{j}))) legend_vec{j}];
            end
            if(j == num_N+1)
                legend_vec{j} = [repmat(' ', 1, 2*max_len-3) legend_vec{j}]; % for infinity
            end
            if(fig_ind < 3) % normalized - give skewness and kurtosis
                legend_vec{j} = [legend_vec{j} ...
                    ' |      ' sprintf('%0.3f', s(j)) '     |    ' sprintf('%0.3f', k(j)) ''];
            else
                legend_vec{j} = [legend_vec{j} ...
                    ' |      ' sprintf('%0.3f', mu(j)) '     |    ' sprintf('%0.3f', sigma(j)^2) ''];
            end
        end
    else
        for j=1:num_N
            legend_vec{j} = ['k=' legend_vec{j}];
        end
    end
    if(fig_ind == 1)
        offset = 1;
    else
        offset = max_len+1;
    end
    if(add_kurtosis)
        if(fig_ind < 3) % normalized
            legend_vec = [[repmat(' ', 1, offset) 'N | Skewness | Kurtosis'] legend_vec]
        else
            legend_vec = [[repmat(' ', 1, offset) 'N | Mean     | Variance'] legend_vec]
        end
    else
%        legend_vec = [[repmat(' ', 1, offset) 'N'] legend_vec]
    end
    legend_handle = legend(legend_vec, 'fontsize', 11, 'fontweight', 'bold', legend_place);
%    set(legend_handle, 'EdgeColor', 'w')
%    set(legend_handle, 'Color', 'none')
    legend boxoff;
    if(legend_place == 1)
        legend_pos = get(legend_handle, 'position'); 
        legend_pos(1:2) = legend_pos(1:2) + [-0.1 -0.1]; 
        set(legend_handle, 'position', legend_pos); 
    end
%    title(['Distribution of QTL values: ' max_flag ' of N Gaussians'], 'fontsize', 14);
    title([repmat(' ', 1, 100) str2title('supp_fig1')], 'fontsize', 16, 'fontweight', 'bold');
    if(fig_ind == 3)
        xlim([-6 6]); % un-normalzied
    else
        xlim([-4 4]);
    end
    ylim([0 max(y_mat(:))+0.01]);
    xlabel('z', 'fontsize', 14, 'fontweight', 'bold');
    ylabel('frequency', 'fontsize', 14, 'fontweight', 'bold');
    if(exist('fig_outfile', 'var'))
        my_saveas(gcf, [fig_outfile fig_str], format_fig_vec);
    end 
end  % loop on figure types



full_figure(0); loglog(1./(2*log(N_vec)), sigma(1:end-1).^2, '.'); hold on;
loglog(1./(2*log(N_vec)), 1./(2*log(N_vec)), 'r'); 
xlabel('2 log(N)'); ylabel('\sigma'); 






% % %
% % %     plot(gumbel_bins_loc, gumbel_hist_normalized, 'k--', 'linewidth', 1);  % Last one: Plot Gumbel distribution
% % %
% % %
% % %
% % % full_figure; % Plot un-normalized trait distributions
% % % for i=1:num_N % loop on number of gaussians
% % %     if(i < 4)
% % %         j=i;
% % %     else
% % %         j=i+1;
% % %     end
% % %     plot(x_vec, y_vec{i}, color_vec(mod_max(j, 96)), 'linewidth', 2); % plot distributions on top of each others. Avoid black
% % % end
% % % legend(legend_vec(1:end-1), 'fontsize', 11);
% % % title(['Distribution of QTL values: Max of N Gaussians (Not Normalized)'], 'fontsize', 14);
% % % xlim([-6 6]);
% % % xlabel('z', 'fontsize', 14, 'fontweight', 'bold');
% % % ylabel('freq.', 'fontsize', 14, 'fontweight', 'bold');
% % % my_saveas(gcf, [fig_outfile '_unnormalized'], format_fig_vec);
% % %
% % % end