% Make dir if it doesn't exist
function my_mkdir(dir_str)
if(~exist(dir_str, 'dir'))
    mkdir(dir_str);
end
