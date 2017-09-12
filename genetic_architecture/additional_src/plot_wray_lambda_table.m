
function plot_wray_lambda_table(lambda_s_vec, lambda_mz_vec, ...
    liab_assum_lambda_s_vec, disease_vec, x_limit)
N = length(lambda_s_vec)

% Flip: lambda_mz on x-axis
loglog(lambda_mz_vec, lambda_s_vec, '*r', 'linewidth', 2); hold on;
loglog(lambda_mz_vec, liab_assum_lambda_s_vec ,'og', 'linewidth', 2);


loglog((1:x_limit).^2, 1:x_limit, ':');
switch  x_limit
    case 100
        x_ticks = [4 5 10 20 30 50 100 200 500 1000 10000];
        y_ticks = [2 3 5 10 15 20 30 50 75 100];
        xlim([4 10000]); ylim([2 100]);
    case 5
        x_ticks = [4 5:5:25];
        y_ticks = [2:5];
        xlim([4 20]); ylim([2 sqrt(20)]);
end


% plot(log10(lambda_s_vec), log10(lambda_mz_vec) ,'*r', 'linewidth', 2);
% plot(log10(lambda_s_vec), log10(liab_assum_lambda_mz_vec) ,'og', 'linewidth', 2);
% plot(log10(1:70), log10((1:70).^2), '--');

legend({'observed \lambda_s', 'expected \lambda_s'}, 2, 'fontsize', 14, 'fontweight', 'bold');
title('Familial risk epidemiological estimates for various diseases (log-scale)');
for i=1:N
    line([lambda_mz_vec(i) lambda_mz_vec(i)], ...
        [lambda_s_vec(i) liab_assum_lambda_s_vec(i)], ...
        'color', 'k', 'linewidth', 2);
    cur_offset = lambda_s_vec(i) / 10;
    new_point = [lambda_mz_vec(i)+10*cur_offset^2 (lambda_s_vec(i) + liab_assum_lambda_s_vec(i))/2 + cur_offset];
    line([lambda_mz_vec(i) new_point(1)], ...
        [(lambda_s_vec(i) + liab_assum_lambda_s_vec(i))/2 new_point(2)], ...
        'color', 'k', 'linewidth', 2, 'linestyle', '--');
    switch i
        case {16,17,18}
            cur_color = 'r';
        otherwise
            cur_color = 'k';
    end
    
    text(new_point(1)+0.01, new_point(2), ...
        disease_vec{i}, 'fontweight', 'bold', 'color', cur_color);
    %    if(i >= N-2)
    %        text(log10(lambda_s_vec(i)), log10(lambda_mz_vec(i)), ...
    %            disease_vec{i}, 'fontsize', 6); % , 'textstyle', 'bold');
    %    end
end

xlabel('\lambda_{MZ} ', 'fontsize', 14, 'fontweight', 'bold');
ylabel('\lambda_s', 'fontsize', 14, 'fontweight', 'bold');

x_labels = num2str_cell(num2cell(x_ticks))
my_xticklabels(x_ticks, x_labels);
y_labels = num2str_cell(num2cell(y_ticks))
my_yticklabels(y_ticks, y_labels, 'YAxisLocation', 'left');


replace_labels = 0;
if(replace_labels) % Replace labels
    x_ticks = get(gca, 'XTick');
    xlabels = num2str_cell(num2cell(x_ticks),2);
    for i=1:length(xlabels)
        xlabels{i} = ['10^{' xlabels{i} '}'];
    end
    my_xticklabels(x_ticks, xlabels);
    y_ticks = get(gca, 'YTick');
    ylabels = num2str_cell(num2cell(y_ticks),2);
    for i=1:length(ylabels)
        ylabels{i} = ['10^{' ylabels{i} '}'];
    end
    my_yticklabels(y_ticks, ylabels, 'YAxisLocation', 'left');
    xlabh = get(gca,'XLabel');
    set(xlabh,'Position',get(xlabh,'Position') - [0 0.1 0])
    ylabv = get(gca,'YLabel');
    set(ylabv,'Position',get(ylabv,'Position') - [0.05 0 0])
end % replace labels



% % % return;
% % %
% % %
% % %
% % %
% % % figure; hold on; plot(lambda_s_vec, lambda_mz_vec ,'*r');
% % % plot(lambda_s_vec, liab_assum_lambda_mz_vec,'+k');
% % % plot(1:70, (1:70).^2);
% % % xlabel('\lambda_s '); ylabel('\lambda_{MZ} ');
% % % legend({'reported', 'expected (LT)', 'expected (mult.)'},1);
% % % title('Epidemiological estimates for various diseases'); %  from Wray''s paper');
% % % ylim([0 2000]);
% % % for i=1:N
% % %     text(lambda_s_vec(i), lambda_mz_vec(i), disease_vec{i});
% % % end
% % % my_saveas(gcf, fullfile(figs_dir, 'lambda_s_vs_lambda_mz_various_diseases'), format_fig_vec);
% % %
