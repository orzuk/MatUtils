% Extract file suffix from a file name 
%
% Input: 
% file_name - the file name (with suffix)
% Output: 
% suffix - just the file suffix
%
function suffix = suffix_from_file_name(file_name)

dot_ind = find(file_name == '.', 1, 'last'); 

if(~isempty(dot_ind))
    suffix = file_name(dot_ind+1:end); 
else % empty: no suffix
    suffix = '';
end

