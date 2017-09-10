% Copy all files in dir1 to dir2 (by default don't overwrite files) 
% 
% The input: 
% dir1 - source directory
% dir2 - destination directory
% overwrite_flag - (optional, default 0) should we overwrite files already present at dir 2
% 
function copy_dir_content(dir1, dir2, overwrite_flag, varargin)

if(~exist('overwrite_flag', 'var'))
    overwrite_flag = 0;
end
if(overwrite_flag)
    copy_files = GetFileNames(dir1);
else
    copy_files = setdiff_dirs(dir1, dir2);
    for i=1:length(copy_files)
        copy_files{i} = fullfile(dir1, copy_files{i}); 
    end
end

for i=1:length(copy_files)
    eval([' copyfile(''' copy_files{i} ''', ''' fullfile(dir2, '.') ''');']);
end

