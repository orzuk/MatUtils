% Take all files appearing in either directory 
function union_files = union_dirs(dir1, dir2, case_sensitive, varargin)
f1 = GetFileNames(dir1);
f2 = GetFileNames(dir2);
if(~exist('case_sensitive', 'var'))
    case_sensitive = 1;
end
if(case_sensitive == 0)
    f1 = upper_cell(f1); f2 = upper_cell(f2); 
end
union_files = union(f1, f2); 



