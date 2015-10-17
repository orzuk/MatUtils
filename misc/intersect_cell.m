% Intersect each element of two cell arrays
function C = intersect_cell(A,B)
n = length(A);
C = cell(n,1);
for i=1:n
    C{i} = intersect(A{i}, B{i});
end


