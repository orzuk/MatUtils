% Match a key string to a value using a cell array of mapping
% The input:
% key - the key to look for the item
% c - cell array of size 2xn. First column is keys and second is values
%
function v = key_to_val(key, c)

if(iscell(key))
    n = length(key);
    v = cell(n,1);
    for i=1:n
        v{i} = key_to_val(key{i}, c);
    end
else
    s = strfind_cell(c(:,1), key);

    if(isempty(s))
        v = [];
        for i=1:size(c,1) % look the other way
            s = strmatch(c(i,1), key);
            if(~isempty(s))
                v = c(i,2);
                break;
            end
        end
    else
        v = c(s,2);
    end
    if(length(v) == 1)
        v = v{1};
    end

end
