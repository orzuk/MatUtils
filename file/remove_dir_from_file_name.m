% Extract only file name from a string containing the file name and its full path
%
% Input: 
% file_name - the file name (full path)
% Output: 
% file_name_no_dir - just the file name (no path)
%
function file_name_no_dir = remove_dir_from_file_name(file_name)

if(iscell(file_name)) % work on many files
    n = length(file_name);
    file_name_no_dir = cell(n,1);
    for i=1:n
        file_name_no_dir{i} = remove_dir_from_file_name(file_name{i});
    end
    return;
end


unix_dir_ind = strfind(file_name, '/');
dir_ind = union(unix_dir_ind, strfind(file_name, '\')); % unite windows and unix

if(~isempty(dir_ind))
    file_name_no_dir = file_name(max(dir_ind)+1:end);
else
    file_name_no_dir = file_name;
end
    
        