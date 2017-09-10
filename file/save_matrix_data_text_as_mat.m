% Save any matrix data in .mat format
%
% Input: 
% input_file - data (with labels) 
% output_file - where to save the data
% empty_numeric_val (optional) what to convert the empty values to (default is zero)
%
function save_matrix_data_text_as_mat(input_file, output_file, empty_numeric_val, varargin)

if(~exist('empty_numeric_val', 'var'))
    empty_numeric_val = 0;
end
R = loadcellfile(input_file); % load expression 

labels = R(1,:); % assume first row is of labels

n = length(labels); % number of labels 
for i=1:n
    labels{i} = replace(labels{i}, ' ', '_'); % get rid of spaces    
    
    full_inds = find(1-isempty_cell(R(2:end,i)));    
    if(min(isnumeric_cell(R(full_inds+1,i)))) % This means that all non-empty cells are numeric
        R(2:end,i) = empty_cell_to_numeric_val(R(2:end,i), empty_numeric_val); % transfer all empty cells to -1
        eval_str = [labels{i} ' = cell2mat(R(2:end,' num2str(i) '));']; % convert to .mat 
    else
        eval_str = [labels{i} ' = R(2:end,' num2str(i) ');']; % leave as string 
    end
    eval(eval_str); 
end

if(~exist('output_file', 'var'))
    output_file = [input_file(1:end-4) '.mat'];
end

save(output_file, 'labels');

for i=1:n
   eval_str = ['save(''' output_file ''',''' labels{i} ''', ''-append'');'];
   eval(eval_str); 
end



