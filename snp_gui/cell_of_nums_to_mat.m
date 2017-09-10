%function mat = cell_of_nums_to_mat(cell_num)
function mat = cell_of_nums_to_mat(cell_num)

% mat = zeros(size(cell_num));
% for i = 1:size(cell_num, 1)
%     for j = 1:size(cell_num, 2)
%         mat(i,j) = str2num(char(cell_num{i,j}));
%     end
% end
% 

f_name = 'LibiTempFile.txt'; % Temporary file
saveCellFile(cell_num, f_name);

mat = load(f_name, 'ascii');

delete(f_name);