% Plot allele frequencies and genetic relative risks for gwas data
function plot_data_stats(data, data_params, output_figs_dir)

AssignGeneralConstants();

binary_inds = find( data.trait_type_num(data_params.unique_inds) == 0);
quant_inds = find( data.trait_type_num(data_params.unique_inds) == 1);

figure; hist(data.RAF, 100); title('Risk-allele Freq. (Controls)'); % Plot risk-allele-frequency histogram
xlabel('Risk Allele Freq.'); ylabel('freq.');
my_saveas(gcf, '../../common_disease_model/figures/database/risk_allele_freq', ...
    format_fig_vec);
do_log = 1;
if(do_log) % perform log transform
    data.OR = log(data.OR); 
    log_label_vec = ' (log-scale)';
else
    log_label_vec = '';
end

figure; hist((data.OR), 100); % Plot odds-ratio histogram
title(['Odds-ratio ' log_label_vec]); 
xlabel(['Odds-ratio' log_label_vec]); ylabel('freq.');
my_loglog('x'); % change x axis ticks
% if(do_log) 
%     xlabels = get(gca, 'XTickLabel');
%     xlabels = num2str(exp(str2num(xlabels)), 3);
%     set(gca, 'XTickLabel', xlabels);
% end
my_saveas(gcf, '../../common_disease_model/figures/database/odds_ratio', ...
    format_fig_vec);

figure; hold on; % plot odds-ratio vs. risk-allele-frequency
data.lethality(data.trait_type_num == 1) = -1; % don't include quantitative traits in lethality yet
for i=1:3
    plot(data.RAF(data.lethality == i-1), data.OR(data.lethality == i-1), ['*' color_vec(i)]); % plot risk-allele-frequency vs. odds-ratio
end
numeric_inds = find(data.trait_type_num == 1);
legend_vec = {'neutral', 'moderate', 'lehtal'};
if(~isempty(numeric_inds))
    plot(data.RAF(numeric_inds), data.OR(numeric_inds), ['.' color_vec(4)]); % plot risk-allele-frequency vs. odds-ratio
    legend_vec = [legend_vec 'QTLs'];
end
label_vec = 'Effect Size (Odds-Ratio)';
if(do_log)
    label_vec = [label_vec ' (log-scale)'];
end
title(['Risk Allele Freq. vs. ' label_vec]);
legend(legend_vec); xlabel('Risk Allele Freq.'); ylabel(label_vec);
if(do_log) % change y axis ticks
    ylabels = get(gca, 'YTickLabel');
    ylabels = num2str(exp(str2num(ylabels)), 3);
    set(gca, 'YTickLabel', ylabels);
end
my_saveas(gcf, '../../common_disease_model/figures/database/odds_ratio_vs_risk_allele_freq', ...
    format_fig_vec);


figure; hold on; % Plot smoothed version of allele odds-ratio (still doesn't work well)
for i=1:3
    [X sort_perm] = sort(data.RAF(data.lethality == i-1));
    Y = data.OR(data.lethality == i-1);
    YY = smooth(X, Y(sort_perm), 'lowess');
    plot(X, YY, color_vec(i));
end
title(['Risk Allele Freq. vs. ' label_vec ' smoothed']);
legend(legend_vec); xlabel('Risk Allele Freq.'); ylabel([label_vec ' (smoothed)']);
my_saveas(gcf, '../../common_disease_model/figures/database/odds_ratio_vs_risk_allele_freq_smoothed', ...
    format_fig_vec);


figure; hold on; % Cumulative plot by RAF
for i=1:3
    X = sort(data.RAF(data.lethality == i-1));
    plot(X, (1:length(X)) ./ length(X), color_vec(i));
end
X = sort(data.RAF(numeric_inds));
plot(X, (1:length(X)) ./ length(X), color_vec(4));
title('Risk Allele Freq. by lethality level');
legend(legend_vec); xlabel('Risk Allele Freq.'); ylabel('cumulative prob.');
my_saveas(gcf, '../../common_disease_model/figures/database/risk_allele_freq_by_lethality', ...
    format_fig_vec);

RAF_bins_starts = [0 0.1 0.3];
RAF_bins_ends = [0.1 0.3 1];
for RAF_ind = 1:length(RAF_bins_starts) % plot in different bins of RAF size
    figure; hold on; % Cumulative plot by OR
    for i=1:3
        bin_inds = find((data.RAF <= RAF_bins_ends(RAF_ind)) & (data.RAF >= RAF_bins_starts(RAF_ind)));
        bin_inds = intersect(bin_inds, find(data.lethality == i-1));
        X = sort(data.OR(bin_inds));
        plot(X, (1:length(X)) ./ length(X),  color_vec(i));
    end
    X = sort(data.OR(numeric_inds));
    plot(X, (1:length(X)) ./ length(X), color_vec(4));
    title(['Risk Allele Freq. by lethality level for SNPs with RAF in [' ...
        num2str(RAF_bins_starts(RAF_ind)) ', ' num2str(RAF_bins_ends(RAF_ind)) ']']);
    legend(legend_vec); xlabel('Odds-Ratio (log)'); ylabel('cumulative prob.');
    my_saveas(gcf, ['../../common_disease_model/figures/database/odds_ratio_by_lethality_RAF_bin_' ...
        num2str(RAF_bins_starts(RAF_ind)) '_' num2str(RAF_bins_ends(RAF_ind))], ...
        format_fig_vec);
end



% Plot histograms showing heritability of diseases
figure; hold on; hist(data_params.h_add ( binary_inds ), 30)
title('Histogram of heritability explained on binary scale');
xlabel('h_{add}'); ylabel('freq.');
my_saveas(gcf, fullfile(output_figs_dir, 'h_additive_binary_scale_hist'), format_fig_vec);

figure; hold on; hist(data_params.h_liab ( binary_inds ), 30)
title('Histogram of heritability explained on liability scale');
xlabel('h_{liab}'); ylabel('freq.');
my_saveas(gcf, fullfile(output_figs_dir, 'h_additive_liability_scale_hist'), format_fig_vec);


figure; hold on; plot(data_params.h_add, data_params.h, '.'); plot(0:0.01:0.5, 0:0.01:0.5, 'r');
xlabel('h_{add}'); ylabel('h_{mult}'); title('Additive vs. multiplicative variance explained estimation');



mu_vec = [0.001 0.01 0.05 0.1 0.2]; % different possible disease prevalences
h_liab_vec = [0.01:0.01:1]; 
figure; hold on; legend_vec = {}; % plot lambda_s vs. heritability: empirical and theoretical
%legend_vec = {}; %  Plot the same for lambda_MZ
for j=1:length(mu_vec) % plot for different mu's
    mu = mu_vec(j);
    [lambda_s_vec_theoretical{j} lambda_mz_vec_theoretical{j}] = ...
        heritability_to_sibling_relative_risk(h_liab_vec, 'liability', mu);
    plot(h_liab_vec, log(lambda_s_vec_theoretical{j}), [color_vec(j) '--']);
    legend_vec{j} = ['\mu=' num2str(mu)]; 
        
end
%title(['Heritability (liability-scale) vs. \lambda_{MZ}']);
data_params.h_familial = data.h_familial(data_params.unique_inds); 
data_params.lambda_s_familial = data.lambda_s_familial(data_params.unique_inds); 
plot_binary_inds = intersect(binary_inds, find(data_params.prevalence ~= 0.01)); % remove all default prevalence
plot_familial_inds = intersect(plot_binary_inds, find(data_params.h_familial >= 0));
plot(data_params.h_liab(plot_binary_inds), log(data_params.lambda_s(plot_binary_inds)), '*'); % Plot lambda_s vs. heritability
plot(data_params.h_familial(plot_familial_inds), log(data_params.lambda_s_familial(plot_familial_inds)), 'r*'); % Plot lambda_s vs. heritability
xlabel('h_{liab}^2'); ylabel('\lambda_s'); title('Heritability and sibling-relative-risk for different diseases');
legend([legend_vec 'diseases (explained)', 'diseases (estimated)' ],2);
for i=1:length(plot_binary_inds) % add text: disease name
   text( data_params.h_liab(plot_binary_inds(i)), log(data_params.lambda_s(plot_binary_inds(i))), ...
       data_params.trait_name{plot_binary_inds(i)} ); 
end
for i=1:length(plot_familial_inds) % add arrow: explained -> estimated 
    arrow([data_params.h_liab(plot_familial_inds(i)), log(data_params.lambda_s(plot_familial_inds(i)))], ...
        [data_params.h_familial(plot_familial_inds(i)), log(data_params.lambda_s_familial(plot_familial_inds(i)))]);
end

% Set Axis to allow more space:
max_x = 1.05 * max( max(data_params.h_liab(plot_binary_inds)), ...
    max(data_params.h_familial(plot_familial_inds)) );
max_y = 1.05 * max( max(log(data_params.lambda_s(plot_binary_inds))), ...
    max(log(max(1,data_params.lambda_s_familial(plot_familial_inds)))) ); % force above 1 to avoid negatives and complex numbers
axis([0 max_x 0 max_y]);

ylabels = get(gca, 'YTickLabel'); % Change labeling at the end
ylabels = num2str(exp(str2num(ylabels)), 3);
set(gca, 'YTickLabel', ylabels);

my_saveas(gcf, fullfile(output_figs_dir, 'heritability_liability_vs_lambda_s_diseases'), format_fig_vec);

