% Sum the elements in each member of a cell array
%
% Input:
% c - cell to sum
% c2 - (optional) sum another cell array
%
% Output:
% M - sum
%
function M = sum_cell(c, c2, dim)

n = length(c);
if(~exist('c2', 'var') || isempty(c2))
    if(~exist('dim', 'var') || isempty(dim))
        M = zeros(n,1);  
        for i=1:n
            M(i) = sum(c{i});
        end
    else
        M = zeros(n,size(c{1}, setdiff(1:2, dim)));  
        for i=1:n
            M(i,:) = sum(c{i}, dim);
        end
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
