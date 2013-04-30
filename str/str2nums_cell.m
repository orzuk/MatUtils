% Convert strings to many-numbers in a cell array
%
% Input:
% s - cell array of strings
% ind - which index to take in each entry (default: take all)
% special_parse_flag - special parsing of percentiles and exponents
%
% Output:
% nc - a cell array of numbers
%
function nc = str2nums_cell(s, ind, special_parse_flag, varargin)
nc = s;
for i=1:size(nc, 1)
    for j=1:size(nc, 2)
        %    do_i = i
        if(isa(nc{i,j}, 'char'))
            if(exist('special_parse_flag', 'var'))
                nc{i,i} = str2nums(s{i,j}, ind, special_parse_flag);
            else
                if(exist('ind', 'var'))
                    nc{i,j} = str2nums(s{i,j}, ind);
                else
                    nc{i,j} = str2nums(s{i,j});
                end
            end
        end
    end
end