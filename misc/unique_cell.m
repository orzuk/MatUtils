% Perform unique to a cell of cell arrays
function [u inds]= unique_cell(c)
n = length(c);
[s lens] = cellfun(@cell2vec, c, 'UniformOutput', false);
[s inds] = unique(s); lens = lens(inds);
u = cellfun(@vec2cell, s, lens, 'UniformOutput', false);




