% Plot the 'enrichment' FDR when comparing two sets,
% where the assumption is that the density of events in both 
% of them is the same
% 
% Input: 
% inds_vec - the instances are assumed to be ordered by rank, where '1'
% means of the set and '0' means of the control
% perm_score_vec - if this is given as input, it is used to sort the inds_vec (top scores are on top)
% alpha - the fraction of the expected instances between the set and the control
% new_fig_flag - flag saying if to start a new figure (default is on)
% plot_curve_flag - flag saying which type of plot to produce (ROC/FDR)
% font_size - plot font size (default is 15) 
% 
% Output: 
% AUC - the area-under-ROC curve (for the case of ROC)
%
function AUC = enrichment_fdr_plot(inds_vec, perm_score_vec, alpha, ...
    legend_vec, new_fig_flag, plot_curve_flag, font_size, varargin)

AssignStatsConstants; % why do we need this? 
AssignGeneralConstants; 
addutils; % just be sure to add utilities 

n = size(inds_vec, 1); % number of different indices
m = size(inds_vec, 2); % number of differnet orientations (what do orientations mean?) 
if( (nargin < 3) || isempty(alpha) )% here alpha not given as input
    alpha = sum(inds_vec) / (n - sum(inds_vec)); % fraction of set from background
end
if(~exist('perm_score_vec', 'var'))
    perm_score_vec = [];
end
if(~isempty(perm_score_vec)) % sort according to the scores vec 
    if(size(perm_score_vec, 2) == 1)
        perm_score_vec = repmat(perm_score_vec, 1, m);
    end
    for i=1:m
        [sorted_score_vec sort_perm] = sort(perm_score_vec(:,i), 'descend');
        inds_vec(:,i) = inds_vec(sort_perm,i); 
    end
end
if(~exist('new_fig_flag', 'var'))
    new_fig_flag = 1;
end
if(new_fig_flag)
    figure;
end
if(~exist('plot_curve_flag', 'var'))
    plot_curve_flag = FDR; 
end
if(~exist('font_size', 'var'))
    font_size = 15;
end
hold on;
switch plot_curve_flag
    case FDR
        title('FDR for given number of true positives', 'fontsize', font_size);
        xlabel('FDR', 'fontsize', font_size); ylabel('log(num. true positives)', 'fontsize', font_size);
        legend_loc = 2; % upper left corner
    case ROC
        title('ROC Curve for given number of motif instances', 'fontsize', font_size);
        xlabel('False Positive Rate', 'fontsize', font_size); ylabel('True Positive Rate', 'fontsize', font_size);
        legend_loc = 4; % lower right corner
end
set(gca, 'fontsize', font_size);
AUC = zeros(m,1); % area under ROC curve variable 
for i=1:m
    set_vec = cumsum(inds_vec(:,i));
    background_vec = cumsum(1-inds_vec(:,i));
    switch plot_curve_flag
        case FDR
            enrichment_fdr_vec = alpha .* background_vec ./ set_vec;
            k = find(set_vec > 0, 1); enrichment_fdr_vec = enrichment_fdr_vec(k:end); set_vec = set_vec(k:end);
            enrichment_fdr_vec = cummin(enrichment_fdr_vec(end:-1:1)); enrichment_fdr_vec = enrichment_fdr_vec(end:-1:1);
            plot(enrichment_fdr_vec, log2(set_vec), color_vec(i), 'linewidth', 3);
        case ROC
            x_vec = background_vec ./ background_vec(end); y_vec = set_vec ./ set_vec(end);
            [x_vec sort_perm] = sort(x_vec); y_vec = y_vec(sort_perm);
            plot(x_vec, y_vec, color_vec(i), 'linewidth', 3);
            AUC(i) = trapz(x_vec, y_vec);            
    end
end
if(exist('legend_vec', 'var'))
    internal_legend_vec_is = legend_vec
    for i=1:m
        legend_vec{i} = [legend_vec{i} ' (AUC ' num2str(AUC(i)) ' )'];
    end
else
    legend_vec = cell(1,m);
    for i=1:m
        legend_vec{i} = ['AUC ' num2str(AUC(i))];
    end
end
legend(legend_vec, legend_loc); % put legend in the best location
   
if(plot_curve_flag == ROC)
    plot(0:0.01:1, 0:0.01:1, 'k.'); % just plot the line y=x
end
    
