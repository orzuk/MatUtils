% Plot nice quantiles
function quantiles_legend_vec = plot_quantiles(x_mat, y_mat, plot_x_vec, log_y_flag)

if(~exist('log_y_flag', 'var') || isempty(log_y_flag))
    log_y_flag = 0;
end

if(size(x_mat,1) == 1)
    x_mat = repmat(x_mat, size(y_mat,1), 1);
end
quantiles_vec = [0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99];
quantiles_legend_vec = num2str(quantiles_vec');
quantile_symbol_vec = {':', ':', '-.', '--', '-'};
quantile_size_vec = [3 3.5 3.5 4 4] .* 0.8;
quantyle_color_vec = {[1 0.7 0.2], [1 0.5 0.15], [1 0.4 0.1], [1 0.15 0.05], 'r'};

z_mat = zeros(size(y_mat,1), length(quantiles_vec)); % plot quantiles 
for j=1:size(y_mat,1) % get quantiles
    for k=1:length(quantiles_vec)
        [~, tmp_quantile_ind] = min(abs(y_mat(j,:) - quantiles_vec(k)));
        z_mat(j,k) = x_mat(j, tmp_quantile_ind);
    end
end

[plot_x_vec sort_perm] = sort(plot_x_vec);  z_mat = z_mat(sort_perm,:); % sort x-vec

for k=1:length(quantiles_vec)
    if(log_y_flag)
        loglog(plot_x_vec, z_mat(:,k), ...
            quantile_symbol_vec{min(k, length(quantiles_vec)+1-k)},  ...
            'color', quantyle_color_vec{min(k, length(quantiles_vec)+1-k)}, ...
            'linewidth', quantile_size_vec( min(k, length(quantiles_vec)+1-k) ));
    else
        semilogx(plot_x_vec, z_mat(:,k), ...
            quantile_symbol_vec{min(k, length(quantiles_vec)+1-k)},  ...
            'color', quantyle_color_vec{min(k, length(quantiles_vec)+1-k)}, ...
            'linewidth', quantile_size_vec( min(k, length(quantiles_vec)+1-k) ));
    end
    hold on;
end
