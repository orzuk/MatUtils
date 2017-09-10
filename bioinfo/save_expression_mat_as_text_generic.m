% Save expression data in .txt format
%
% Input:
% labels - a row cell array
% cell_array - the first few columns
% data - the data (in numerical format)
% output_file - name of file to output
% in_matlab_flag - if to do stuff in matlab (slower but numbers look nicer)
% output_format - (new, optional) also enable different formats such as .gct
%
function save_expression_mat_as_text_generic(labels, cell_array, data, ...
    output_file, in_matlab_flag, output_format, varargin)

n = size(data,1); % num rows
m = size(data,2);  % num data columns
k = size(cell_array, 2); % num extra columns

if(~exist('in_matlab_flag', 'var') || isempty(in_matlab_flag))
    in_matlab_flag = 0;
end
if(~exist('output_format', 'var'))
    output_format = suffix_from_file_name(output_file); 
end
if(in_matlab_flag) % do everything in matlab - could be slow for large files
    R = cell(n+1,m+k);
    R(1,m+k-length(labels)+1:end) = labels; % missing labels are assumed to be at the beginning
    if(~isempty(cell_array))
        R(2:end,1:k) = cell_array;
    end
    for i=1:n
        for j=1:m
            R{i+1,j+k} = num2str(data(i,j)); % make string - looks better when saving 
        end
    end
    
    switch output_format
        case 'gct'
            if(k == 1) % duplicate first column
                R = [R(:,1) R]; % duplicate first line
            end
            R = [cell(2,size(R,2))' R']';
            R{1,1} = '#1.2';
            R{2,1} = n; R{2,2} = m;
            R{3,1} = 'Name'; R{3,2} = 'Description'; 
    end        
    savecellfile(R, output_file);
    
else  % use some scripting - should be faster for large files
    %    data = [(1:n)' data];
    %    cell_array = [mat2cell(1:n, 1, ones(n,1))' cell_array];

%     save_str_flag = 1; % nicer saving, but larger running time (for big matrices)
%     if(save_str_flag)
%         data = num2str(data);
%     end
    save([output_file(1:end-4) '.data.txt'], 'data', '-ASCII'); % save the data
    %     for i=1:n
    %         for j=1:m
    %             R{i,j} = data(i,j);
    %         end
    %     end
    %     savecellfile(R, [output_file(1:end-4) '.data.txt']);
    
    if(~isempty(cell_array))
        cell_array = replace_cell(cell_array(:,1), '', '-');
        savecellfile(cell_array, [output_file(1:end-4) '.all_labels.txt']); % save the labels
        system(['python join_files_columns.py ' [output_file(1:end-4) '.all_labels.txt '] ...
            [output_file(1:end-4) '.data.txt'] ' ' output_file ]);     % combine both files
    else
        system(['cp ' output_file(1:end-4) '.data.txt ' output_file ]);
    end
    
    savecellfile(labels, [output_file(1:end-4) '.all_labels.txt'], ' ');
    system(['cat ' [output_file(1:end-4) '.all_labels.txt '] ' ' output_file ' > ' [output_file(1:end-4) '.data.txt'] ]);
    system(['python change_file_separators.py ' [output_file(1:end-4) '.data.txt '] output_file ]);
    my_delete([output_file(1:end-4) '.all_labels.txt ']);
    my_delete([output_file(1:end-4) '.data.txt ']);
end





