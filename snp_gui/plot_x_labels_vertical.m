%function plot_x_labels_vertical(x_loc_vec, y_loc, labels, font_size)
function plot_x_labels_vertical(x_loc_vec, y_loc, labels, font_size)

if(nargin < 4)
    font_size = 7;
end
set(gca,'xtick',[]);

text(x_loc_vec, repmat(y_loc, size(x_loc_vec)), labels,...
    'rotation', 270, 'FontSize', font_size);
