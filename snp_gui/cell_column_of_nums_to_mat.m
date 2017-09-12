%function mat = cell_column_of_nums_to_mat(cell_num)
function mat = cell_column_of_nums_to_mat(cell_num)

mat = [];
if(size(cell_num,2)==1)
    t=char(cell_num);
    mat = str2num(t);
end