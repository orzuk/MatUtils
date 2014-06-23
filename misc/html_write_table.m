% Write a table in html format
function R_tab = html_write_table(R, header_flag, sortable_flag)

% % % /* Sortable tables */
% % % table.sortable thead {
% % %     background-color:#eee;
% % %     color:#666666;
% % %     font-weight: bold;
% % %     cursor: default;
% % % }

[num_rows num_columns] = size(R); 

if(~exist('header_flag', 'var') || isempty(header_flag))
    header_flag = 0; 
end
R_tab{1} = '<table border="1"  class="sortable">';
for i=1:num_rows 
    R_tab{i+1} = '<tr>';    
    if(header_flag && (i==1)) % make header boldface
        for j=1:num_columns
            R_tab{i+1} =  [R_tab{i+1} '<th bgcolor="#CCCCCC">' R{i,j} '</th> ']; % light gray
        end
    else
        for j=1:num_columns
            R_tab{i+1} =  [R_tab{i+1} '<td>' R{i,j} '</td> '];
        end
    end
    R_tab{i+1} = [R_tab{i+1} ' </tr>'];
end
R_tab{end+1} = '</table>'; 

% <tr>
% <td>row 1, cell 1</td>
% <td>row 1, cell 2</td>
% </tr>
% <tr>
% <td>row 2, cell 1</td>
% <td>row 2, cell 2</td>
% </tr>
% </table>
