%function new_cell_str = add_pref_to_cell_str(cell_str, str)
function new_cell_str = add_pref_to_cell_str(cell_str, str)

if(size(cell_str, 1)==1) cell_str = cell_str'; end

num_entries = size(cell_str, 1);
cell_to_add = cell(num_entries, 1);
cell_to_add = assign_str_to_indices_in_cell(cell_to_add, str, [1:num_entries]);
new_cell_str = cellstr([char(cell_to_add) char(cell_str)]);

