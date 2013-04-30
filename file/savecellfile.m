% Save a cell array into a tab-delimited .txt file (Written by Assif Yitzhaky)
%
% Input:
% table - the cell array to be saved
% outfile - where to save the cell
% delimiter - what delimiter to use when saving to file (default is tab '\t')
% save_mode - whether to save strings without the '%s' (default, 0) or with
% (1). In the default state, '\t' will be converted into tabs, and '\n' to new lines
%
% New addition: some cells in the input are typically empty. We now convert
% them to empty strings to avoid crashing problems. 
% 
function savecellfile(table, outFile, delimiter, save_mode, varargin)

table = empty_cell_to_empty_str(table); % convert empty cells to empty strings 

if((nargin < 3) || isempty(delimiter))
    delimiter = '\t';
end
if(nargin < 4) 
    save_mode = 0;
end
if(~isnumeric(save_mode))
    save_mode= 0;
end
fout=fopen(outFile, 'wt');

% saving_lines = size(table,1)
for i=1:size(table, 1)  % loop on lines 
   % i_is = i
    if(mod(i, 1000) == 0)
        sprintf('save line %ld out of %ld',  i, size(table,1))
    end
    for j=1:size(table, 2)-1
        if isempty(table{i,j})
            fprintf(fout, delimiter);
            continue;
        elseif isnumeric(table{i,j})
            specifier=['%d' delimiter];
        else
            specifier=['%s' delimiter];
        end
        fprintf(fout, specifier, table{i,j});
    end
    if isnumeric(table{i, size(table, 2)}) % last column treated differently. Why? 
        specifier='%d';
    else
        specifier='%s';
    end
    if(~isempty(strmatch('%s', specifier))  && (~save_mode) && isempty(strfind(table{i, size(table, 2)}, 'href=' ))) % we don't use the %s - to avoid tabs being lost .. - but then the problem is that \n are turned into new lines .. and we also lose '%' at the end!
        fprintf(fout, table{i, size(table, 2)} );
    else
        fprintf(fout, specifier, table{i, size(table, 2)} );
    end
    fprintf(fout, '\n');
end
fclose(fout);
