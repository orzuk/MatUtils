function move_additional_src()

path(path, fullfile('..', 'src'));

func_list = depfun('SNP_GUI');

last_ind = strfind(func_list, 'src\');

ctr=1;
for i=1:length(func_list)
    if(~isempty(last_ind{i}))
        src_func_list{ctr} = func_list{i}(last_ind{i}+4:end);
        ctr=ctr+1;
    end
end

s = dir(fullfile('..', 'src', '*.m'));
files_in_src_folder = {};
ctr=1;
for i = 1:length(s)
    files_in_src_folder{ctr} = s(i).name;
    ctr=ctr+1;
end

files_to_move = setdiff(files_in_src_folder, src_func_list);

for i = 1:length(files_to_move)
    i
   mv_str = ['move "' fullfile('..', 'src', files_to_move{i}) '" ' '"' files_to_move{i} '"']
   system(mv_str);
end