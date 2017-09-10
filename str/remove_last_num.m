% Remove last number from a string
function s_no_last_num = remove_last_num(s)

if(iscell(s))
    s_no_last_num = cellfun(@remove_last_num, s, 'uniformoutput', false);
else
    num_inds = ( (s <= '9') & (s >= '0'));
    n = length(s);
    for j=n:-1:1
        if(~num_inds(j))
            break;
        end
    end
    s_no_last_num = s(1:j);
end
