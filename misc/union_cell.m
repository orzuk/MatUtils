% Unite each element of two cell arrays
%
% Input:
% A - first cell array
% B - second cell array
%
% Output:
% C - union cell array
%
function C = union_cell(A,B, union_str, varargin)
n = size(A,1); m = size(A,2);
if(~exist('union_str', 'var') || isempty(union_str))
    union_str = '';
end
if(exist('B', 'var'))
    C = cell(n,m);
    for i=1:n
        for j=1:m
            if(isempty(union_str) || (~isempty(A{i,j}) && ~isempty(B{i,j})))
                C{i,j} = union(A{i,j}, B{i,j}, union_str);
            else
                if(isempty(A{i,j}))
                    C{i,j} = B{i,j};
                else
                    C{i,j} = A{i,j};
                end
            end
        end
    end
else % perform union of just A
    C = A{1,1};
    for i=1:n
        for j=1:m
            if(isempty(union_str) || (~isempty(A{i,j}) && ~isempty(C)))
                C = union(C, A{i,j}, union_str);
            else
                if(isempty(C))
                    C = A{i,j};
                end
            end
        end
    end
    %    C = unique(cell2vec(A));
end



