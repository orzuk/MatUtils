% Check for each element of a cell array if it is numberic
function M = isnumeric_cell(c)

n = length(c);
M = zeros(n,1); 
for i=1:n
    M(i) = isnumeric(c{i});
end
