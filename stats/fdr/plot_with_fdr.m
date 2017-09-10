% Plot a set of p-values with the FDR line
function plot_with_fdr(pvals_vec, fdr_alpha,  title_str, new_fig_flag, varargin)

if(exist('new_fig_flag', 'var'))
    figure;
end
hold on; plot(sort(pvals_vec), '.');
n = length(pvals_vec);
x_vec = [1:n];
plot(x_vec, x_vec .* fdr_alpha ./ n, 'r');

[num_rejected, fdr_vec, idx] = fdr(pvals_vec, fdr_alpha);
num_rejected_is = num_rejected
if(~exist('title_str', 'var'))
    title_str = [];
end
title([title_str ' rej. ' num2str(num_rejected) ' out of ' num2str(n) ' hypot. at FDR = ' num2str(fdr_alpha)]);



