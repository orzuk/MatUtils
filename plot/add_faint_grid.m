% Add a faint grid to figure
% Input: 
% gray_level - how dark is grid (higher means darker?)
% keep_axes - keep axis original color (black) 
% 
function add_faint_grid(gray_level, keep_axes)

if(~exist('gray_level', 'var') || isempty(gray_level))
    gray_level = 0.7;
end
if(~exist('keep_axes', 'var') || isempty(keep_axes))
    keep_axes = 1;
end
grid on;
if(keep_axes)
Caxes = copyobj(gca, gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); % keep axis colors intact 
end

set(gca,'Xcolor',[gray_level, gray_level, gray_level]);
set(gca,'Ycolor',[gray_level, gray_level, gray_level]);
set(gca,'GridLineStyle',':');

