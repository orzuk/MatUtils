% Make a string nice for plotting
%
% Input:
% str_in - original sring
%
% Output:
% str - output string (prettier, no underscores, greek laters)
%
function str = str2title(str_in)

if(iscell(str_in))
    str = cellfun(@str2title, str_in, 'uniformoutput', false);
else
    str_in = [' ' str_in];
    str_in = strrep(str_in, '_', '-');
    greek = {'alpha', 'beta', 'gamma', 'delta', 'chi', 'xi', 'pi', 'mu', ...
        'sigma', 'omega', 'lambda', 'epsilon', 'phi', 'theta'};
    [~, J] = strfind_cell(lower(str_in), greek);
    
    if(~isempty(J))
        J_vec = max(1,J-1);
        rep_inds = find(str_in(max(1,J-1)) ~= '\');
        if(~isempty(rep_inds) && (max(J_vec(rep_inds)) > 1))
            rep_inds = setdiff(rep_inds, ...
                find((str_in(max(1,J-1)) >= 'a') & (str_in(max(1,J-1)) <= 'z')));
            rep_inds = setdiff(rep_inds, ...
                find((str_in(max(1,J-1)) >= 'A') & (str_in(max(1,J-1)) <= 'Z')));
        end
        
        J = sort(J(rep_inds)); % f = find(str_in(inds-1) ~= '\');
        if(~isempty(J))
            n = length(str_in); J_comp = ones(length(J)+n,1); J_comp(J+(0:length(J)-1)) = 0;
            str = repmat('\', 1, n+length(J)); str(J_comp) = str_in;
        else
            str = str_in;
        end
    else
        str = str_in;
    end
    str = strrep(str, '\-', '_'); % replace back to keep original under-scores marked with '\_' 
    str = str(2:end); 
end
