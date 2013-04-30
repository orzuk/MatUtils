% Checks the premissions of files in a certain directory (works only with unix).
% In the directory field one can give also a description of the files - e.g.
% '/data/*.txt' and the function will look only for such files
%
% Input: 
% file_dir_str - string with the file directory, e.g. '/data/*.txt'
% with_dir - flag saying to return also the directory (default is false)
% 
% Output:
% file_names - names of files in dir
% file_permissions - matrix representing permission of all files
%
function [file_names file_permissions] = GetFilePermissions(files_dir_str, with_dir, varargin)

if(~exist('with_dir', 'var') || isempty(with_dir))
    with_dir = 0;
end

file_names = GetFileNames(files_dir_str, with_dir); 
n = length(file_names); 
file_permissions = zeros(n,9); % 9 permissions: read/write/execute X user/group/everybody
system(['ls -l ' files_dir_str ' > tmp_ls_res.txt']);


R = loadcellfile('tmp_ls_res.txt', [], 9);
R = strsplit_cell(R(2:end), ' ');
for i=1:length(R)
    j = strmatch(R{i}{end}, remove_dir_from_file_name(file_names), 'exact');
    if(~isempty(j))
        file_permissions(j,:) = ~(R{i}{1}(2:10) == '-');
    end
end
system('rm tmp_ls_res.txt'); % remove the tmp file 



