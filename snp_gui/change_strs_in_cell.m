%function out_cell = change_strs_in_cell(in_cell, str_to_change_cell, new_str)
function out_cell = change_strs_in_cell(in_cell, str_to_change_cell, new_str)

out_cell = in_cell;
num_str_to_change  = size(str_to_change_cell, 1);

all_str_ind = [];
for i = 1:num_str_to_change 
    str = char(str_to_change_cell(i));
    str_ind = strcmp(in_cell, str);
    str_ind = find(str_ind==1);
    all_str_ind = [all_str_ind str_ind'];
end

if(length(all_str_ind)>0)
    out_cell(all_str_ind) = {new_str};
end
