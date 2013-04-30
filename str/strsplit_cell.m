% Split each string in a cell array 
function parts = strsplit_cell(c, delim)
n = length(c); 
parts = cell(n,1);
for i=1:n
    parts{i} = strsplit(c{i}, delim); 
end


