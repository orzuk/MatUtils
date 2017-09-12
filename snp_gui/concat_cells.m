%function new_cell = concat_cells(cell1, cell2, dim)
function new_cell = concat_cells(cell1, cell2, dim)

new_cell = [];
if(nargin==2)
    if(size(cell2,2) == size(cell1,2))
        new_cell = cell(size(cell1,1)+size(cell2,1), size(cell2,2));
        new_cell(1:size(cell1,1),:) = cell1;
        new_cell(size(cell1,1)+1:end,:) = cell2;
    elseif(size(cell2,1) == size(cell1,1))
        new_cell = cell(size(cell2,1), size(cell1,2)+size(cell2,2));
        new_cell(:, 1:size(cell1,2)) = cell1;
        new_cell(:, size(cell1,2)+1:end) = cell2;
    end
elseif(nargin==3)
    if(dim==2)
        new_cell = cell(size(cell1,1)+size(cell2,1), max(size(cell2,2), size(cell1,2)));
        new_cell(1:size(cell1,1),1:size(cell1,2)) = cell1;
        new_cell(size(cell1,1)+1:end,1:size(cell2,2)) = cell2;
    elseif(dim==1)
        new_cell = cell(max(size(cell2,1), size(cell1,1)), size(cell1,2)+size(cell2,2));
        new_cell(1:size(cell1,1), 1:size(cell1,2)) = cell1;
        new_cell(1:size(cell2,1), size(cell1,2)+1:end) = cell2;
    end
end
