% Plot comparison of selection coefficients fitted for different populations and different genes
%
% Input:
% FittedParams- structure representing selection parameters
% exome_data_figs_dir - where to save output figures
%
function plot_fitted_selection_parameters(FittedParams, exome_struct, exome_data_figs_dir)

if(ischar(FittedParams)) % load from file
    FittedParams = load(FittedParams);
end
% get population names
num_populations = min(length(exome_struct.populations), 7);

S_corr_mat = corr(FittedParams.s_MLE_vec); % compute pairwise correlation between populations
% Filter unfitted values
for i=1:num_populations
    for j=1:num_populations
        good_inds{i,j} = find((max(FittedParams.s_MLE_vec(:,[i j])') < 1) & (min(FittedParams.s_MLE_vec(:,[i j])') > 0));
        S_corr_mat(i,j) = corr(FittedParams.s_MLE_vec(good_inds{i,j},i), FittedParams.s_MLE_vec(good_inds{i,j},j));
    end
end

% Plot figure: heatmap
figure; imagesc_with_labels(S_corr_mat, exome_struct.populations, exome_struct.populations);
%colorbar;
%set(gca, 'xticklabels', exome_struct.populations);
%set(gca, 'yticklabels', exome_struct.populations);
my_saveas(gcf, fullfile(exome_data_figs_dir, 's_population_correlation'), {'epsc', 'pdf', 'jpg'});

figure; % Make scatter-plot
[ha, pos] = tight_subplot(num_populations, num_populations, [0.04 .04],[0.04 0.04],[0.04 0.04]);
for i=1:num_populations
    for j=(i+1):num_populations % plot only above diagonal
        axes(ha(num_populations*(i-1)+j));
        %        subplot(num_populations, num_populations, (i-1)*num_populations+j);
        %        good_inds = find(max(abs(FittedParams.s_MLE_vec(:,i)), abs(FittedParams.s_MLE_vec(:,j)))>0);
        loglog(FittedParams.s_MLE_vec(good_inds{i,j},i), FittedParams.s_MLE_vec(good_inds{i,j},j), '.');
        if(j == 1) % population name
            ylabel(exome_struct.populations{i}, 'fontsize', 7);
        end
        if(i == 1) % population name
            title(exome_struct.populations{j}, 'fontsize', 7);
        end
        set(gca, 'FontSize', 7);
    end % loop on j
    % on diagonal we plot histograms of each population
    axes(ha(num_populations*(i-1)+i));
    %    subplot(num_populations, num_populations, (i-1)*num_populations+i);
    hist(log(FittedParams.s_MLE_vec(unique(cell2mat(good_inds(i,:))),i)), 100); % use log-scale?
    if(i == 1) % population name
        title(exome_struct.populations{i}, 'fontsize', 7);
    end
    set(gca,'xtick',[]); set(gca,'ytick',[]);
        set(gca, 'FontSize', 7);
    % below diagonal we can plot something else:
end
set(0, 'defaultFigureRenderer', 'painters'); % set render to show dots
my_saveas(gcf, fullfile(exome_data_figs_dir, 's_population_scatter'), {'epsc', 'pdf', 'jpg'});

% Figure: Compare our fitted selection with the fitted coefficients of Samocha et al.




