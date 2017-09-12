%function ret = is_file_pres(f_name)
function ret = is_file_pres(f_name)

fid = fopen(f_name);
if(fid ~= -1)
    ret = 1;
else
    ret = 0;
end

if(ret)
    fclose(fid);
end