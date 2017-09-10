% remove strings from each member of a cell-array of strings
% 
% Input: 
% c - cell-array of strings 
% s - string to remove 
% 
% Output: 
% c - new string without s
% 
function c = strdiff_cell(c, s)
for i=1:length(c)
    c{i} = strdiff(c{i}, s); 
end

