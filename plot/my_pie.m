% Prepare a pie figure but include only positive values
function my_pie(x, labels, varargin)

pos_inds = find(x > 0);
if(exist('labels', 'var'))
else
    labels = num2str_cell(num2cell(round(1000 * x(pos_inds) ./ sum(x(pos_inds)))./10));
    labels = cellfun(@strcat, labels, repcell('%', length(labels), 1), 'uniformoutput', false);
end
pie(x(pos_inds), labels);
    