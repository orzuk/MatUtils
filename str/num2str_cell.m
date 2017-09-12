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
    if(~exist('precision', 'var') || isempty(precision))
        for i=1:size(nc,1)
            for j=1:size(nc,2)
                if(isnumeric(nc{i,j}))
                    s{i,j} = num2str_delim(nc{i,j}, delimiter);
                end
            end
        end
    else
        for i=1:size(nc,1)
%            i_is = i 
%            size_is = size(nc,1)
            for j=1:size(nc,2)
                if(isnumeric(nc{i,j}))
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




