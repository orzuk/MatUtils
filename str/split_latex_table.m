% Split a table to a few table (so it doesn't go out of pages)
%
% Input:
% latex_tab - table (in cell-array of strings)
% tab_header - before the table (including caption)
% tab_footer - at the end of the table (including caption)
% num_rows_in_page - max. # of rows to keep in one page before splitting
% split_format - -1 ??? counter subtracting to keep table number constant
%
% Output:
% latex_tab_split - vector of table split into pages
%
function latex_tab_split = split_latex_table(latex_tab, ...
    tab_header, tab_footer, num_rows_in_page, split_format)

if(~exist('split_format', 'var') || isempty(split_format)) % counter setting
    split_format = [];
end
num_rows = size(latex_tab, 1)-1; % remove first line 
num_tables = ceil((num_rows-1) / num_rows_in_page);

latex_tab_split = [];
for i=1:num_tables
    cur_rows = [1 ((i-1)*num_rows_in_page+2):min(i*num_rows_in_page+1, num_rows)];
    
    latex_block_to_add = ...
        [vec2row(tab_header) vec2row(latex_tab(cur_rows)) vec2row(tab_footer)];
    for j=1:length(latex_block_to_add) % make sure labels are different
        latex_block_to_add{j} = strrep(latex_block_to_add{j}, '\label{', ...
            ['\label{split_' num2str(i) '_']);
    end
    
    if(~isempty(split_format))
        if(i > 1)
            latex_block_to_add = [['\addtocounter{table}{' num2str(split_format) '}'] latex_block_to_add]; % add counter update
            label_ind = strfind_cell(latex_block_to_add, '\label');
            if(~isempty(label_ind))
                label_ind_in_line = strfind(latex_block_to_add{label_ind}, '\label')
                label_ind_in_line2 = strfind(latex_block_to_add{label_ind}, '}')
                label_ind_in_line2 = label_ind_in_line2(find(label_ind_in_line2>label_ind_in_line, 1))
                latex_block_to_add{label_ind} = ...
                    [latex_block_to_add{label_ind}(1:label_ind_in_line-1) ...
                    latex_block_to_add{label_ind}((label_ind_in_line2+1):end)]; % remove label
            end
            caption_ind = strfind_cell(latex_block_to_add, '\caption');
            if(~isempty(caption_ind))
                latex_block_to_add{caption_ind}  = [latex_block_to_add{caption_ind}(1:end-1) ' (cont.) }']; % add '(cont.)'
                square_parent_ind = strfind(latex_block_to_add{caption_ind}, '[')% remove first part of caption []
                square_parent_ind2 = strfind(latex_block_to_add{caption_ind}, ']')% remove first part of caption []
                if(~isempty(square_parent_ind2))
                    latex_block_to_add{caption_ind} = [latex_block_to_add{caption_ind}(1:square_parent_ind) ...
                        latex_block_to_add{caption_ind}((square_parent_ind2):end)];
                end
            end
        end
    end
    latex_tab_split = [vec2row(latex_tab_split) latex_block_to_add];
end
latex_tab_split = vec2column(latex_tab_split);





