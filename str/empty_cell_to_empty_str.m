% Convert empty cells in a cell array to empty strings
function c_str = empty_cell_to_empty_str(c)

c_str = c;
for i=1:size(c_str,1)
    for j=1:size(c_str,2)
        if(isempty(c_str{i,j}))
            c_str{i,j} = '';
        end
    end
end

    
