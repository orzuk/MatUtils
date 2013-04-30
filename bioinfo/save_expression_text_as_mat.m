% Save expression data in .mat format
%
% Input:
% input_file - data file with expression matrix
% output_file - where to save the data
% overwrite_flag - whether to over write an existing .mat file (default is yes)
%
function save_expression_text_as_mat(input_file, output_file, overwrite_flag, varargin)

if(~exist('output_file', 'var') || isempty(output_file))
    output_file = file_name_to_mat(input_file);
end
if(~exist('overwrite_flag', 'var'))
    overwrite_flag = 1;
end

if( (~exist(output_file, 'file')) || overwrite_flag )
    R = loadcellfile(input_file); % load expression
    first_row = min(strfind_cell(R(:,1), '_'));
    first_column = find(isnumeric_cell(R(4,:)), 1);

    num_genes = size(R,1)-first_row+1; % number of genes
    num_samples = size(R,2)-first_column+1;
    labels = R(first_row-1,first_column:end);
    affy_id = R(first_row:end,1); % first column is affy ids
    gene_symbol = R(first_row:end,2); % second column is gene symbol
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

    save(output_file, 'data', 'labels', 'affy_id', 'gene_symbol', 'num_genes', 'num_samples');
end
