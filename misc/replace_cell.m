% Replace every occurance of a in x by b for a cell x
%
% The input: 
% x - cell array
% a - value to be replaced
% b - value to replace
% 
% The output: 
% x - new cell array after replacing
% 
function x = replace_cell(x,a,b)

if(isa(a, 'char')) % replace a char
    v = vec2row(my_strmatch(a, x, 'exact'));
    if(~isempty(v))
        for i= v
            x{i} = b;
        end
    end
else
    v = vec2row(strfind_cell(x, a));
    if(~isempty(v))
        for i= v
            x{i} = b;
        end
    end
end
