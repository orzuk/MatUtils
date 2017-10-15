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
num_populations = min(length(exome_struct.populations), 4); 
S_corr_mat = corr(FittedParams.s_MLE_vec); % compute pairwise correlation between populations 

% Plot figure: heatmap 
figure; imagesc_with_labels(S_corr_mat, exome_struct.populations, exome_struct.populations); 
%colorbar;
%set(gca, 'xticklabels', exome_struct.populations);  
%set(gca, 'yticklabels', exome_struct.populations);  
my_saveas(gcf, fullfile(exome_data_figs_dir, 's_population_correlation'), {'epsc', 'pdf', 'jpg'}); 

figure; % Make scatter-plot
for i=1:num_populations
    for j=(i+1):num_populations % plot only above diagonal
        subplot(num_populations, num_populations, (i-1)*num_populations+j);
        good_inds = find(max(abs(FittedParams.s_MLE_vec(:,i)), abs(FittedParams.s_MLE_vec(:,j)))>0);
        plot(FittedParams.s_MLE_vec(good_inds,i), FittedParams.s_MLE_vec(good_inds,j), '.');
        if(j == 1) % population name
            ylabel(exome_struct.populations{i}, 'fontsize', 14);
        end
        if(i == num_populations) % population name
            xlabel(exome_struct.populations{j}, 'fontsize', 14);
        end
    end
    % on diagonal we plot histograms of each population
    subplot(num_populations, num_populations, (i-1)*num_populations+i);
    hist(log(abs(nonzeros(FittedParams.s_MLE_vec(:,i)))), 100); % use log-scale? 
    
    % below diagonal we can plot something else:
end
set(0, 'defaultFigureRenderer', 'painters'); % set render to show dots 
my_saveas(gcf, fullfile(exome_data_figs_dir, 's_population_scatter'), {'epsc', 'pdf', 'jpg'}); 

% Figure: Compare our fitted selection with the fitted coefficients of Samocha et al.




