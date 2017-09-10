%function new_str_cell = erase_backslash_cell(str_cell)
function new_str_cell = erase_backslash_cell(str_cell)

num_raws = size(str_cell, 1);
flag = 0;
if(num_raws == 1)
    str_cell = str_cell';
    num_raws = size(str_cell, 1);
    flag = 1;
end
new_str_cell = str_cell;
for i = 1:num_raws
    str = char(str_cell(i,:));
    new_str = str;
    backslash_ind = findstr('_', new_str);
    new_str(backslash_ind) = '-';
    new_str_cell{i,1} = new_str;
end

if(flag)
    new_str_cell = new_str_cell';
end