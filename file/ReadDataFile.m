% Read a file in a generic format. 1st line gives field names, other lines give the data
%
% Input:
% data_file - name of data (.txt) file
% output_file - where to save .mat file
% cell_to_mat - flag saying if to convert all numeric fields to matrices (default is on)
% skip_lines - skip first lines of file
% delimiter - how to separate fields 
%
% Output:
% S - struct with fields read
% R - 'raw' data (sometimes needed too)
% R_skipped - 'raw' data skipped
%
function [S R R_skipped] = ReadDataFile(data_file, output_file, cell_to_mat, skip_lines, delimiter, varargin)

%loading_file = data_file
if(~exist('cell_to_mat', 'var') || isempty(cell_to_mat)) % default is converting to .mat 
    cell_to_mat = 1;
end
if(~exist('skip_lines', 'var') || isempty(skip_lines))
    skip_lines = 0;
end
if(~exist('delimiter', 'var') || isempty(delimiter))
    delimiter = [];
end
if(is_mat_file(data_file))
    load(data_file);
else
    to_mat = (cell_to_mat >= 0); % need to put -1 if we want it not to work
    R = loadcellfile(data_file, to_mat, delimiter);
    if(exist('output_file', 'var')) % save in .mat format 
        if(isempty(output_file))
            output_file = file_name_to_mat(data_file);
        end
        save([remove_suffix_from_file_name(output_file) '_raw.mat'], 'R');
    end
end


% Allow skip_lines to be a character. We avoid the lines containing this character
if(ischar(skip_lines))
    skip_lines_inds = strfind_cell(R(:,1), skip_lines); 
    keep_lines_inds = setdiff(1:length(R(:,1)), skip_lines_inds);
    skip_lines = keep_lines_inds(1)-1;    
end
    
% field_names = strrep_cell(strrep_cell(R(1,:), ' ', '_'), '.', '_'); %
% avoid spaces and dots
R_skipped = R(1:skip_lines,:); % NEW! Save also skipped indices
R = R((1+skip_lines):end,:);
if(isempty(R)) % empty file
    S = []; 
    return;
end
field_names = strrep_cell(R(1,:), {' ', '.', ',', '-', '(', ')', '{', '}', '[', ']', ...
    '\', '/', '#', '%', ':', '+', '?', '"', '<', '>', '=', char(160), char(194)}, '_');% avoid spaces,dots, parenthesis and slashes

S = [];
for i=1:length(field_names)
    run_field = i
    field_name_is = field_names{i}
    if(isempty(field_names{i}))
        continue;
    end
    if( (field_names{i}(1) == '_') || ...
            ((field_names{i}(1) <= '9') && (field_names{i}(1) >= '0')) )
        field_names{i} = ['XXX' field_names{i}];
    end
    eval_str = ['S.' field_names{i} ' = R(2:end,' num2str(i) ');']
    eval(eval_str);
end

if(nargout == 1) 
    clear R; % save memory (unless we want this as output) 
end
if(cell_to_mat > 0) % convert to num fields with only numeric values
    for i=1:length(field_names)
        if(isempty(field_names{i})) % don't read empty fields
            continue;
        end
        eval_str = ['num_flag = all(isnumeric_cell_local(S.' field_names{i} ...
            ') | isempty_cell_local(S.' field_names{i} '));']; eval(eval_str);
        if(num_flag)
            eval_str = ['S.' field_names{i} ...
                ' = cell2mat(empty_cell_to_numeric_val_local(S.' field_names{i} ', 0));']; eval(eval_str);
        end
    end
end
if(exist('output_file', 'var'))
    if(isempty(output_file))
        output_file = file_name_to_mat(data_file);
    end
    save(output_file, '-struct', 'S');
end

% Check for each element of a cell array if it is numberic
function M = isnumeric_cell_local(c)

n = length(c);
M = zeros(n,1); 
for i=1:n
    M(i) = isnumeric(c{i});
end

% Convert empty cells or strings in a cell array to a value (default is zero)
function c_str = empty_cell_to_numeric_val_local(c, v, varargin)

if(~exist('v', 'var'))
    v = 0;
end

c_str = c;
for i=1:length(c_str)
    if(isempty(c_str{i}))
        c_str{i} = v;
    end
end

% Check if cell entries are empty 
function v = isempty_cell_local(c)
n = length(c); v = zeros(n,1);
for i=1:n
   v(i) = isempty(c{i}); 
end
