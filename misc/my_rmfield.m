% Remove field but check that it exists first
function y = my_rmfield(x, s)

f = find(isfield(x, s));

if(~isempty(f)) % some fields to remove
    if(iscell(s)) % take only relevant fields 
        s = s(f);
    end
    y = rmfield(x, s);
else
    y = x;
end

