% Extracts only the directory from a string containing the file
% name and its full path
% Input: 
% file_name - the file name (full path)
% Output: 
% file_dir - just the directory without the file name 
%
function file_dir = dir_from_file_name(file_name)

unix_dir_ind = strfind(file_name, '/');
dir_ind = union( union(unix_dir_ind, strfind(file_name, '|')), strfind(file_name, '\') ); % unite windows and unix

if(~isempty(dir_ind))
    file_dir = file_name(1:max(dir_ind));
else
    file_dir = ''; % no directory specified
end
    
