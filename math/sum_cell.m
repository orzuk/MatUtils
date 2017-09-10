% Sum the elements in each member of a cell array
function M = sum_cell(c, c2)

n = length(c);
if(nargin == 1)
    M = zeros(n,1);
    for i=1:n
        M(i) = sum(c{i});
    end
else % here add two cells
    M = c;
    if(iscell(c2)) % add two cell arrays
        for i=1:n
            M{i} = M{i} + c2{i};
        end
    else % add same scalar/vector each time
        for i=1:n
            M{i} = M{i} + c2;
        end
    end
end
