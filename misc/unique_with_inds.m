% Perform unique but also keep indices. Just call 'get_duplicates' 
function [vals, inds, num_dups] = unique_with_inds(v)

[vals, inds, num_dups] = get_duplicates(v);
