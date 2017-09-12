% Get a file with the same file name as input file but different suffix
function new_file_name = file_name_to_new_suffix(file_name, suffix)
new_file_name = [remove_suffix_from_file_name(file_name) '.' suffix];



