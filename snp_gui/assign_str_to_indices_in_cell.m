%function assigned_cell = assign_str_to_indices_in_cell(cell_array, str, ind_vec)
function assigned_cell = assign_str_to_indices_in_cell(cell_array, str, ind_vec)

assigned_cell = cell_array;
assigned_cell(ind_vec) = {str};