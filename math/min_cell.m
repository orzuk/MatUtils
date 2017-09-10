% Get the minimum element in each member of a cell array
function [m, I] = min_cell(c)

n = length(c);
m = zeros(n,1); I=m;
for i=1:n
    [m(i), I(i)] = min(c{i}(:));
end
