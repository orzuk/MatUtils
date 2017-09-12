% Upper for a cell, deals with empty cells (from Tal Shay)
function A = my_upper(cellstr)

I = find(~cellfun('isempty', cellstr));
B = upper(cellstr(I));
A = cell(size(cellstr));
A(I) = B;
