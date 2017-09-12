% Unite columns from different files to the same file
% 
% Input: 
% file_names - list of files
% col_vec - which columns to take from each file 
% output_file - name of output file
% labels_vec - ?? 
% 
% Output: 
% T - ? 
% 
function T = merge_file_columns(file_names, col_vec, output_file, labels_vec)

if(~iscell(file_names))
    file_names = GetFileNames(file_names, 1);
end
n = length(file_names);
if(length(col_vec) == 1)
    col_vec = repmat(col_vec, n, 1);
end
if(~iscell(col_vec))
    col_vec = num2cell(col_vec);
end
T = {}; ctr = 0; use_labels_vec = {};
for i=1:n
    if(exist(file_names{i}, 'file'))
        R = loadcellfile(file_names{i});
        T = [T R(:,col_vec{i})];
        if(exist('labels_vec', 'var')) % update labels_vec
            labels_is = labels_vec
            ccc = col_vec{i}+ctr
            use_labels_is = use_labels_vec
            use_labels_vec = [use_labels_vec vec2row(labels_vec((ctr+1) : (ctr+length(col_vec{i})) ))]
            ctr = ctr + length(col_vec{i})
        end
    end
end
if(exist('labels_vec', 'var')) % add labels on top row
    T = [vec2column(use_labels_vec) T']';
end
savecellfile(T, output_file);

