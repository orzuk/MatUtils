% Convert empty cells or strings in a cell array to a value (default is zero)
function c_str = empty_cell_to_numeric_val(c, v, varargin)

if(~exist('v', 'var'))
    v = 0;
end

c_str = c;
for i=1:length(c_str)
    if(isempty(c_str{i}))
        c_str{i} = v;
    end
end
