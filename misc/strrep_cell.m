% Replace all occurances of a string a in a cell array c by b
% 
% Input: 
% c - cell array
% a - string to be replaced
% b - replacing string
% 
% Output: 
% c - new cell array after replaced strings
% 
function c = strrep_cell(c, a, b)
if(iscell(a)) % replace multiple strings 
    for i=1:length(a)
        c = strrep_cell(c,a{i},b);
    end
else
    for i=1:size(c,1) % enable two-dimensional cell array 
        for j=1:size(c,2)
            if(ischar(c{i,j})) % replace only in strings (not numbers) 
                c{i,j} = strrep(c{i,j}, a, b);
            end
        end
    end
end
