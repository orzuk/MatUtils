% Plot a cumulative sum of a vector
% 
% The input: 
% x_vec - data vector
% q - where to put quantile line
% new_fig_flag - whether to open a new figure
% title_str - what to put in title (cell array of 3 strings)
% xlabel_str - what to put in x axis
% ylabel_str - what to put in y axis
% outfile - where to save results
% fig_format_vec - what format to save results in
%
function Dummy = plot_cumsum(x_vec, q, new_fig_flag, title_str, xlabel_str, ylabel_str, outfile, fig_format_vec, varargin)

Dummy = []; 
cumsum_vec = cumsum(sort(x_vec, 'descend'));
if(~exist('new_fig_flag', 'var')) % default is to generate a new figure 
    new_fig_flag = 1;
end
if(new_fig_flag)
    figure; hold on; 
end

plot(cumsum_vec, 'linewidth', 2); 
quant_val = cumsum_vec(end) * q;
[dummy quant_ind] = min(abs(cumsum_vec - quant_val));
line([ quant_ind quant_ind], [0 quant_val], 'linewidth', 2, 'linestyle', '--', 'color', 'r');
line([ 0 quant_ind], [quant_val quant_val], 'linewidth', 2, 'linestyle', '--', 'color', 'r');

if(exist('title_str', 'var'))
    title_str = [title_str{1} ' - ' num2str(quant_ind) ' ' title_str{2} ' take ' num2str(q) ' of ' title_str{3}]; 
    my_title(title_str);
end
if(exist('xlabel_str', 'var'))
    my_xlabel(xlabel_str);
end
if(exist('ylabel_str', 'var'))
    my_ylabel(ylabel_str);
end

if(exist('outfile', 'var'))
    if(~iscell(fig_format_vec))
        tmp = fig_format_vec; fig_format_vec = {tmp};
    end
    for i=1:length(fig_format_vec)
        saveas(gcf, outfile, fig_format_vec{i});
    end
end

    
    

