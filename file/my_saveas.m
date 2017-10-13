% Like saveas but make directory if it's non-existance
%
% Input:
% x - variable to save (e.g. gcf)
% file_name - where to save it
% fomrat - what format to put the fig in (can be cell array)
%
function my_saveas(x, file_name, format)

if(iscell(format))
    for i=1:length(format)
        my_saveas(x, file_name, format{i});
    end
else
    my_mkdir(dir_from_file_name([file_name '.']));  % why mkdir here?
    saving_figure_file = file_name
    saving_format = format
    if(strmatch('eps', format))
        format_str = 'eps';
    else
        format_str = format;
    end
    saveas(x, [file_name '.' format_str], format);
end

