% Get a .mat file with the same file name as input file
function mat_file_name = file_name_to_mat(file_name)
mat_file_name = [remove_suffix_from_file_name(file_name) '.mat'];



