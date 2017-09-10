% Split each string in a cell array
% Input:
% c - cell array with strings
% delim - delimiter for spring split
%
% Output:
% parts - cell array with splitted strings
%
function parts = strsplit_cell(c, delim)
[m, n] = size(c);
parts = cell(m, n);
for i=1:m
    for j=1:n
        parts{i,j} = strsplit(c{i,j}, delim);
    end
end


