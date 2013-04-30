%function ret_cell = assign_empty_str_to_empty_cells(input_cell)
function ret_cell = assign_empty_str_to_empty_cells(input_cell)

ret_cell = input_cell;
m = size(input_cell, 1);
n = size(input_cell, 2);
for i = 1:m
    for j = 1:n
        if(isempty(input_cell{i,j}))
            ret_cell{i,j} = '';
        end
    end
end
