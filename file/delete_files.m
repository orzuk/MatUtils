% Delete all files matching a string in sub-directories
function delete_files(root_dir, file_str, depth)

if(~exist('depth', 'var'))
    depth = 999;
end
g = get_subdirectories(root_dir, depth); 

for i=1:length(g)
   delete(fullfile(g{i}, file_str)); % try all at once
   f = GetFileNames(fullfile(g{i},  file_str), 1); % try separately
   for j=1:length(f)
       delete(f{j});
   end
end
