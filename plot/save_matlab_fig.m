% Save a matlab .fig file in other formats 
function save_matlab_fig(matlab_fig_file_str)

matlab_fig_file_names = GetFileNames(matlab_fig_file_str, 1);
for i=1:length(matlab_fig_file_names)
    open(matlab_fig_file_names{i});
    my_saveas(gcf, remove_suffix_from_file_name(matlab_fig_file_names{i}), ...
        {'epsc', 'pdf', 'jpg'}); % save same name with different extension
    close all;
end
