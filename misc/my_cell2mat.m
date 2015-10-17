% Convert cell to mat (check first if already .mat)
function c = my_cell2mat(x)
if(iscell(x))
    c = cell2mat(x);
else
    c = x;
end
