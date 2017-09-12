% Save expression data in .mat format
%
% Input:
% input_file - data file with expression matrix
% output_file - where to save the data
% overwrite_flag - whether to over write an existing .mat file (default is yes)
% num_cols_labels - number of labels for the columns (e.g. genes)
% num_rows_labels - number of labels (for rows (e.g. samples)
%
% Output: 
% data - matrix with numeric expression values
% labels - label of each sample
% genes - names of all genes
% num_genes - how many genes
% num_samples - how many samples
%
function [data labels genes num_genes num_samples] = ...
    save_expression_text_as_mat_generic(input_file, output_file, overwrite_flag, ...
     num_cols_labels, num_rows_labels, varargin)

if(~exist('output_file', 'var') || isempty(output_file))
    output_file = file_name_to_mat(input_file);
end
if(~exist('overwrite_flag', 'var') || isempty(overwrite_flag))
    overwrite_flag = 1;
end
if(~exist('num_cols_labels', 'var') || isempty(num_cols_labels))
    num_cols_labels = 1;
end
if(~exist('num_rows_labels', 'var') || isempty(num_rows_labels))
    num_rows_labels = 1;
end


if( (~exist(output_file, 'file')) || overwrite_flag )
    R = loadcellfile(input_file); % load expression
    first_row = num_rows_labels+1; % min(strfind_cell(R(:,1), '_'));
    first_column = num_cols_labels+1; % find(isnumeric_cell(R(4,:)), 1);

    num_genes = size(R,1)-first_row+1; % number of genes
    num_samples = size(R,2)-first_column+1;
    labels = R(first_row-1,first_column:end);
    genes = R(first_row:end,1); % first column is affy ids
    for i=1:size(R,1) % get rid of all empty stuff
        for j=1:size(R,2)
            if(isempty(R{i,j}))
                R{i,j} = -1;
            end
        end
    end
    data = single(cell2mat(R(first_row:end,first_column:end))); % get data

    if(~exist('output_file', 'var'))
        output_file = [input_file(1:end-4) '.mat'];
    end

    save(output_file, 'data', 'labels', 'genes', 'num_genes', 'num_samples');
end
