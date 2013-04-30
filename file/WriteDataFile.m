% Write a file in a generic format. 1st line gives field names, other lines give the data
%
% Input:
% data_struct - name of data structure
% output_file - where to save .txt file
% cell_to_mat - flag saying if to convert all numeric fields to matrices (default is on)
% skip_lines - skip first lines of file
%
% Output:
% R - 'raw' data (sometimes needed too)
%
function R = WriteDataFile(data_struct, output_file, cell_to_mat, skip_lines, varargin)

%loading_file = data_file
if(~exist('cell_to_mat', 'var'))
    cell_to_mat = 1;
end
if(~exist('skip_lines', 'var') || isempty(skip_lines))
    skip_lines = 0;
end

% field_names = strrep_cell(strrep_cell(R(1,:), ' ', '_'), '.', '_'); %
% avoid spaces and dots

%%% line skipping R = R((1+skip_lines):end,:);


%field_names = strrep_cell(R(1,:), {' ', '.', '-', '(', ')', '{', '}', '[', ']', ...
%    '\', '/', '#', '%', char(160), char(194)}, '_');% avoid spaces,dots, parenthesis and slashes

field_names = fieldnames(data_struct);



for i=1:length(field_names)
    %    eval_str = ['S.' field_names{i} ' = R(2:end,' num2str(i) ');'];
    if(i == 1) % here set R
        eval_str = ['num_elements = length(data_struct.' field_names{i} ');'];
        eval(eval_str);
        R = cell(num_elements+1, length(field_names));
    end
    R{1,i} = field_names{i};
    eval_str = ['R(2:end,' num2str(i) ') = data_struct.' field_names{i} ';'];
    
    to_cell_str = ['if(isnumeric(data_struct.' field_names{i} ')) data_struct.' field_names{i} ...
        ' = num2str_cell(num2cell(data_struct.' field_names{i} ')); end'];
    eval(to_cell_str);
    eval(eval_str);
end

if(cell_to_mat > 0)
    for i=1:length(field_names)
        eval_str = ['num_flag = all(isnumeric_cell_local(data_struct.' field_names{i} '));']; eval(eval_str);
        if(num_flag)
            eval_str = ['data_struct.' field_names{i} ...
                ' = cell2mat(empty_cell_to_numeric_val_local(data_struct.' field_names{i} ', 0));']; eval(eval_str);
        end
    end
end
if(exist('output_file', 'var'))
    if(isempty(output_file))
        output_file = file_name_to_mat(data_file);
    end
    savecellfile(R, output_file, [], 1); % save to output file. Default is with last line with specifier
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
