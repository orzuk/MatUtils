% Check if a file is .mat file
function ret = is_mat_file(file_name)
switch suffix_from_file_name(file_name)
    case 'mat'
        ret = 1;
    otherwise
        ret = 0;
end
