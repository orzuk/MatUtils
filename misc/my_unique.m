% Unique of a cell, deals with empty cells (from Tal Shay)
function U = my_unique(L)

nonemptyL = L(find(~cellfun('isempty', L)));
U = unique(nonemptyL);