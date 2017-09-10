% Plot the relation between heritability on liability scale
% and sibling-relative risk under the assumption of a multiplicative model
% (is this valid?).
% The goal is to see if heritability and relative risk estimates from
% different papers are consistent.
%
% Input: 
% mu_vec - disease prevalence / trait's means
% lambda_s_vec - sibling relative risks
% h_liab_vec - heritability (on liability scale) 
% labels_vec - labels (disease names)
% fig_outfile - where to save output
% 
function plot_heritability_parameters(mu_vec, ...
    lambda_s_vec, h_liab_vec, labels_vec, fig_outfile)

AssignGeneralConstants();

% mu_vec = [0.001 0.01 0.05 0.1 0.2];
% mu_vec = [0.004]; % T1D prevalence

tol = 0.02; % takes time ..
all_h_liab_vec = tol:tol:1; full_figure; % hold on; % Plot lambda_s vs. heritability (on liability scale)
all_lambda_s_vec = {}; all_lambda_mz_vec = {};
for i=1:length(mu_vec)
    [all_lambda_s_vec{i} all_lambda_mz_vec{i}] = heritability_to_sibling_relative_risk(all_h_liab_vec, 'liability', mu_vec(i));
    plot(all_h_liab_vec, log(all_lambda_s_vec{i}), color_vec(i));
%    plot(all_h_liab_vec, log(sqrt(all_lambda_mz_vec{i})), [color_vec(i) '--']); % No need for sqrt. This is wrong
    plot(h_liab_vec(i), log(lambda_s_vec(i)), [color_vec(i) '*']); % plot a point estimate for a disease
    text(h_liab_vec(i), log(lambda_s_vec(i)), labels_vec{i}, 'color', color_vec(i));
end
% my_loglog('y');
y_ticks = get(gca, 'YTick');
if(length(y_ticks) >= 9)
    y_ticks = interp(y_ticks, 2); % refine resolution
end
set(gca, 'YTick', y_ticks);
y_labels = num2str(exp(y_ticks'), 3);
%ylabels = get(gca, 'YTickLabel');
%ylabels = num2str(exp(str2num(ylabels)), 3);
set(gca, 'YTickLabel', y_labels);

tmp_legend_vec = num2str_cell(num2cell(mu_vec*100)); legend_vec = cell(2*length(mu_vec),1);
for i=1:length(tmp_legend_vec)
    legend_vec{2*i-1} = [labels_vec{i} ' \mu = ' tmp_legend_vec{i} '%, h^2 = ' ...
        num2str(100*h_liab_vec(i)) '%, \lambda_s=' num2str(lambda_s_vec(i))];
    legend_vec{2*i} = '';
end
legend(legend_vec,2, 'fontsize', 9); xlabel('h^2 (liability)'); ylabel('\lambda_s');
title('Heritability on liability scale vs. sibling relative risk for different prevalence levels (assuming liability-threshold model)');

if(exist('fig_outfile', 'var'))
    my_saveas(gcf, fig_outfile, format_fig_vec);
end

figure; hold on; % Plot lambda_mz vs. lambda_s (under liability threshold)
plot(log(sqrt(all_lambda_mz_vec{1})), log(all_lambda_mz_vec{1}), 'k--');
plot(log((all_lambda_mz_vec{1}+1)./2), log(all_lambda_mz_vec{1}), 'k:');
for i=1:length(mu_vec)
    plot(log(all_lambda_s_vec{i}), log(all_lambda_mz_vec{i}), color_vec(i));
    %    loglog((all_lambda_s_vec{i}), (all_lambda_mz_vec{i}), color_vec(i));
end

axis(log([1 10 1 100]))
legend(['Risch Mult.' 'Risch Add.' legend_vec(1:2:end)']', 4); 
xlabel('\lambda_s'); ylabel('\lambda_{MZ}');
%my_loglog('xy')
y_ticks = get(gca, 'YTick');
if(length(y_ticks) > 8)
    y_ticks = interp(y_ticks, 2); % refine resolution
else
    y_ticks = sort([y_ticks y_ticks(1:end-1) + diff(y_ticks) ./ 2]);
end
set(gca, 'YTick', y_ticks);
y_labels = num2str(exp(y_ticks'), 3);
set(gca, 'YTickLabel', y_labels);
x_ticks = get(gca, 'YTick');
if(length(x_ticks) > 8)
    x_ticks = interp(x_ticks, 2); % refine resolution
else
    x_ticks = sort([x_ticks x_ticks(1:end-1) + diff(x_ticks) ./ 2]);
end
set(gca, 'XTick', x_ticks);
x_labels = num2str(exp(x_ticks'), 3);
set(gca, 'XTickLabel', x_labels);
if(exist('fig_outfile', 'var'))
    my_saveas(gcf, [fig_outfile '_lambda_s_vs_lambda_mz'], format_fig_vec);
end


figure; hold on; % Plot another form of isolines: prevalence vs. lambda_s (under liability threshold)
all_prevalence_vec = tol:tol:0.5;
h_vec = [0.01 0.05 0.1 0.2 0.5 0.8];
for i=1:length(h_vec)
    all_lambda_s_vec = heritability_to_sibling_relative_risk(h_vec(i), 'liability', all_prevalence_vec);
    plot(all_prevalence_vec, log(all_lambda_s_vec), color_vec(i));
end
y_ticks = get(gca, 'YTick');
if(length(y_ticks) > 8)
    y_ticks = interp(y_ticks, 2); % refine resolution
end
set(gca, 'YTick', y_ticks);
y_labels = num2str(exp(y_ticks'), 3);
set(gca, 'YTickLabel', y_labels);

tmp_legend_vec = num2str_cell(num2cell(h_vec)); legend_vec = cell(length(h_vec),1);
for i=1:length(tmp_legend_vec)
    legend_vec{i} = [' h^2 = ' tmp_legend_vec{i}];
end
title('Heritability on liability scale for different sibling relative risk and prevalence levels (assuming liability-threshold model)');
legend(legend_vec); xlabel('\mu (prevalence)'); ylabel('\lambda_s');
if(exist('fig_outfile', 'var'))
    my_saveas(gcf, [fig_outfile '_prevalence_vs_lambda_s'], format_fig_vec);
end


