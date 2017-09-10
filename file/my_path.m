% Add directory to path, but check first that it exists (to avoid matlab warnings)
function my_path(path_dir)

if(exist(path_dir, 'dir'))
    path(path, path_dir);
end
