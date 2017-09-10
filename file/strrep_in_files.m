% Replace all occurances of a string in a set of file by another string
% Be careful while using! the function steps on the input files (might ruin them)
% No spaces can be allowed in file names
%
% Input:
% files_dir - directory with input files
% s1 - string to be replaced
% s2 - string to replace it
% output_files_dir - output directory (default - same as input. File is modified) enables saving results in a different file
%
function strrep_in_files(files_dir, s1, s2, output_files_dir, varargin)
if(~exist('output_files_dir', 'var') || isempty(output_files_dir))
    output_files_dir = files_dir;
end

if(iscell(s1)) % replace many strings
    strrep_in_files(files_dir, s1{1}, s2{1}, output_files_dir); % first time copy to a new directory 
    for i=2:length(s1) % next times work on copied version 
        replace_word_ind = i
        replace_word = s1{i}
        strrep_in_files(fullfile(output_files_dir, remove_dir_from_file_name(files_dir)), ...
            s1{i}, s2{i}, output_files_dir);
    end
else
    file_names = GetFileNames(files_dir, 1);
    num_files = length(file_names);
    
    for i=1:num_files
        output_file_name = fullfile(dir_from_file_name(output_files_dir), remove_dir_from_file_name(file_names{i}));
        sed_str = ['sed ''s/' s1 '/' s2 '/g'' ' file_names{i} ' > ' output_file_name '.tmp'];
        system(sed_str);
        copyfile([output_file_name '.tmp'], output_file_name);
        delete([output_file_name '.tmp']);
    end
end

