% Prepare a two-tail quantile figure
% 
% Input: 
% fig_output_file - where to save figure 
% t - quantiles for each tail 
% n - number of individuals in each tail
% r - number of rare alleles in each tail
% beta - effect size 
%
function make_QTL_two_tails_figure(fig_output_file, t, n, r, beta)

AssignGeneralConstants;
full_figure;

x_vec = -5:0.01:5; y_vec = normpdf(x_vec); z_vec = normpdf(x_vec - beta); % Draw stadard Gaussian
s = norminv(t); % get quanitle thresholds
left_inds = find(x_vec < s(1)); right_inds = find(x_vec > s(2));

plot(x_vec, y_vec, 'linewidth', 3); % plot gaussian
plot(x_vec, z_vec, 'linewidth', 3, 'linestyle', '-', 'color', 'r'); % plot shifted gaussian


line([s(1) s(1)], [0 max(normpdf(s(1)), normpdf(s(1)-beta))], ...
    'color', 'k', 'linewidth', 3, 'linestyle', '--');
line([s(2) s(2)], [0 max(normpdf(s(2)), normpdf(s(2)-beta))], ...
    'color', 'k', 'linewidth', 3, 'linestyle', '--');
text(s(1), -0.01, 's_1', 'fontsize', 14, 'fontweight', 'b'); 
text(s(2), -0.01, 's_2', 'fontsize', 14, 'fontweight', 'b'); 

%jbfill(x_vec(right_inds), y_vec(right_inds), zeros(1, length(right_inds)), 'm'); 
jbfill(x_vec(left_inds), z_vec(left_inds), zeros(1, length(left_inds)), 'r'); 
jbfill(x_vec(right_inds), z_vec(right_inds), zeros(1, length(right_inds)), 'r'); 
%jbfill(x_vec(left_inds), zeros(1, length(left_inds)), z_vec(left_inds),  'm'); 
jbfill(x_vec(left_inds), y_vec(left_inds), zeros(1, length(left_inds))); 
jbfill(x_vec(right_inds), y_vec(right_inds), zeros(1, length(right_inds))); 


legend({'controls', 'carriers'}, 'fontsize', 14); 
xlabel('z', 'fontsize', 14); ylabel('density', 'fontsize', 14); 

arrow([2 0.05], [2+beta 0.05]); %arrow([beta/3 0.05], [0 0.05]);
text(2+beta/3, 0.06, '\beta',  'fontsize', 18, 'fontweight', 'b');
text(-3, 0.05, 'r_1/n_1', 'fontsize', 18, 'fontweight', 'b'); 
text(3, 0.05, 'r_2/n_2', 'fontsize', 18, 'fontweight', 'b'); 

if(~isempty(fig_output_file))
    my_saveas(gcf, fig_output_file, format_fig_vec);
end

