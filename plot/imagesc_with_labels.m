% Draw an imagesc but change the axes labels
%
% Input:
% data - data matrix
% x_vec - vector of labels for x axis
% y_vec - vector of labels for y axis
% x_N - height of x labels (default is zero)
% y_N - width of y labels (default is zero) 
%
function imagesc_with_labels(data, x_vec, y_vec, x_N, y_N, varargin)

if(~exist('x_N', 'var') || isempty(x_N))
    x_N = 0;
end
if(~exist('y_N', 'var') || isempty(y_N))
    y_N = 0;
end

h = imagesc_nan(data); colorbar; % draw picture

n = size(data,1);
set(h, 'ydata', [n:-1:1]); % flip the y axis

x_tics = get(gca, 'xtick');
if(length(x_vec) == size(data,2))
    x_tics = 1:size(data,2);
    set(gca, 'xtick', x_tics);
end
x_labels = num2str_cell(vec2column(x_vec(x_tics)));

horz_flag = 0;
if(horz_flag) % horizonal labels
    set(gca,'xtickLabel',x_labels)
else
    for i=1:length(x_labels)
        text(i, x_N, x_labels{i}, 'rotation', 90, 'linewidth', 2);
    end
end

if(exist('y_vec', 'var') && ~isempty(y_vec))
    y_tics = get(gca, 'ytick');
    if(length(y_vec) == size(data,1))
        y_tics = 1:size(data,1);
        set(gca, 'ytick', y_tics);
    end
    y_labels = num2str_cell(vec2column(y_vec(y_tics(end:-1:1))));
%     for i=1:length(y_labels)
%         text(y_N, i, y_labels{i}, 'linewidth', 2);
%     end
    set(gca,'ytickLabel',y_labels);
end


 