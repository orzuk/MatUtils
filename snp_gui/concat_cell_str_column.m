%function cell_str = concat_cell_str_column(cell_str1, cell_str2)
function cell_str = concat_cell_str_column(cell_str1, cell_str2)

%cell_str = cellstr([char(cell_str1) char(cell_str2)]);
cell_str = strcat(cell_str1, cell_str2);

if(size(cell_str1,1) ~= size(cell_str,1) | size(cell_str1,2) ~= size(cell_str,2))
    cell_str = cell_str';
end