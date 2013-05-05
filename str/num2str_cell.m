% Convert numbers to string in a cell array
%
% Input:
% nc - numeric cell array
% precision - rounding
% delimiter - how to seperate different numbers
%
% Output:
% s - string cell array
%
function s = num2str_cell(nc, precision, delimiter, varargin)
% if(~exist('precision', 'var') || isempty(precision))
%     precision = -1;
% end
if(~exist('delimiter', 'var') || isempty(delimiter))
    delimiter = '';
end
if(iscell(nc)) % input is cell
    s = nc;
    for i=1:size(nc,1)
        for j=1:size(nc,2)
            if(isnumeric(nc{i,j}))
                if(~exist('precision', 'var') || isempty(precision))
                    s{i,j} = num2str_delim(nc{i,j}, delimiter);
                else
                    s{i,j} = num2str_delim(nc{i,j}, delimiter, precision);
                end
            end
        end
    end
else % input is vector
    if(~exist('precision', 'var') || isempty(precision))
        s = strtrim(mat2cell(num2str(vec2column(nc)), ones(1, length(nc))));
    else
        s = strtrim(mat2cell(num2str(vec2column(nc), precision), ones(1, length(nc))));
    end
    if(isrow(nc))
        s=s';
    end
end




