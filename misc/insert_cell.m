% Insert the cell array d in the i'th place at the cell array c.
% If d is only one element it can be a non-cell
function c = insert_cell(c, d, i)

if(~iscell(d))
    d = {d};
end
if(isrow(c))
    c = [c(1:i-1) d c(i:end)];
else
    c = [c(1:i-1)' d' c(i:end)']';
end
