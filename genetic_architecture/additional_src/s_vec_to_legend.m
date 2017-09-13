% Get nice legend vectors 
function s_legend_vec = s_vec_to_legend(s_vec)

s_legend_vec = [repmat('s= 10^{', length(s_vec), 1) num2str(log10(abs(s_vec')),3) ...
    repmat('}', length(s_vec), 1)];
s_legend_vec = cellstr(s_legend_vec);
s_legend_vec = strrep_cell(s_legend_vec, ' ', '');
s_legend_vec{s_vec == 0} = 's= 0'; % fix s=0
