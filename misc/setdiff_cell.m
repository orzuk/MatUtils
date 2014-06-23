% Setdiff each element of cell arrays
function C = setdiff_cell(A,B)
n = length(A);
C = cell(n,1);
for i=1:n
    C{i} = setdiff(A{i}, B{i});
end

