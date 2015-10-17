% Convert a cell with numbers to a matrix 
function mat = cell_of_nums_to_mat(cell_num)

f_name = 'temp_file_stam.txt';
saveCellFile(cell_num, f_name);

mat = load(f_name, 'ascii');
delete(f_name);