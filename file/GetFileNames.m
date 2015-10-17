% Get the file names matching a certain regular expression in a directory.
% In the directory field one can give also a description of the files - e.g.
% '/data/*.txt' and the function will look only for such files
%
% Input:
% file_dir_str - string with the file directory, e.g. '/data/*.txt'
% with_dir - flag saying to return also the directory (default is false)
% min_date - New! get only files newer than this date
%
% Output:
% file_names - names of files in dir
% file_sizes - sizes of files
%
function [file_names file_sizes] = GetFileNames(files_dir_str, with_dir, min_date, varargin)

if(~exist('with_dir', 'var') || isempty(with_dir))
    with_dir = 0;
end

f = dir(files_dir_str);

if(exist('min_date', 'var')) % filter by date
    num_files = length(f); new_inds = zeros(num_files,1);
    for i=1:num_files
        new_inds(i) = (etime( datevec(f(i).date), datevec(min_date) ) >= 0);
    end
    f = f(find(new_inds));
end

num_files = length(f); file_names = {}; file_sizes = zeros(1,num_files);
bad_files = []; % files without permission 
file_names = cell(1,num_files); 
for i=1:num_files
    if(with_dir)
        file_names{i} = fullfile(dir_from_file_name(files_dir_str),  f(i).name);
    else
        file_names{i} = f(i).name;
    end
    if(isempty(f(i).bytes)) % this means there is no permission!
        bad_files = [bad_files i];
    else
        file_sizes(i) = f(i).bytes;
    end
end
if(~isempty(bad_files))
    good_files = setdiff(1:num_files, bad_files);
    file_names = file_names(good_files); file_sizes = file_sizes(good_files);
end
