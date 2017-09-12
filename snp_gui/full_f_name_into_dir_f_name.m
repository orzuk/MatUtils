%function [dir, f_name] = full_f_name_into_dir_f_name(full_f_name)
function [dir, f_name] = full_f_name_into_dir_f_name(full_f_name)

ind = strfind(full_f_name, '\');
dir = full_f_name(1:max(ind));
f_name = full_f_name(max(ind)+1:end);
