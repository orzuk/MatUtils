% Return the indices in a cell-array
% of strings in which a given query string s is found
%
% Input:
% c - the cell array
% s - the string to look for
%
% Output:
% I - the indices in the cell array where we find the string s
% J - second indices in the cell array where we find the string s (for 2dim cells)
% K - the place in each word where the string was (first) found
%
function [I, J, K] = strfind_cell(c, s)

I = []; J = []; K = [];
if(iscell(c)) % c is cell array. Can be two dimensional
    for i=1:size(c,1)
        for j=1:size(c,2)
            k = strfind(c{i,j}, s);
            if(~isempty(k))
                I = [I i]; J = [J j]; K = [K k];
            end
        end
    end
    if(my_width(c) == 1) % one-dim vector
        I = max(I,J); J = ones(1,length(I));
    end
else % here s is cell array
    for i=1:length(s)
        j = strfind(c, s{i});
        if(~isempty(j))
            I = [I i]; J = [J j];
        end
    end
end


