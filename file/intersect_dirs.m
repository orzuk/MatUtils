% Check which files are shared between two directories
function inter_files = intersect_dirs(dir1, dir2, case_sensitive, varargin)
f1 = GetFileNames(dir1);
f2 = GetFileNames(dir2);

if(~exist('case_sensitive', 'var'))
    case_sensitive = 1;
end
if(case_sensitive == 0)
    f1 = upper_cell(f1); f2 = upper_cell(f2); 
end
inter_files = intersect(f1, f2); 



