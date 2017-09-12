% Extract file name from a string containing file name and its suffix
%
% Input: 
% file_name - the file name (with suffix)
% Output: 
% file_name_no_suffix - just the file name (no suffix)
%
function file_name_no_suffix = remove_suffix_from_file_name(file_name)

dot_ind = find(file_name == '.', 1, 'last'); 

if(~isempty(dot_ind))
    file_name_no_suffix = file_name(1:dot_ind-1); 
else
    file_name_no_suffix = file_name;
end

