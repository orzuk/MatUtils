%function new_cell_str = add_suff_to_cell_str(cell_str, str)
function new_cell_str = add_suff_to_cell_str(cell_str, str)

num_entries = size(cell_str, 1);
new_cell_str = cell_str;
for i = 1:num_entries
    temp = char(cell_str(i));
    new_cell_str{i} = [temp str];
end
