% Delete file if it exists
function my_delete(file_name)
if(exist(file_name, 'file'))
    delete(file_name);
end
