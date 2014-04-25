% Add a faint grid to figure
% Input: 
% gray_level - how dark is grid (higher means darker?)
% 
function add_faint_grid(gray_level)

if(~exist('gray_level', 'var') || isempty(gray_level))
    gray_level = 0.7;
end
grid on;

Caxes = copyobj(gca,gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); % keep axis colors intact 

set(gca,'Xcolor',[gray_level gray_level gray_level]);
set(gca,'Ycolor',[gray_level gray_level gray_level]);
set(gca,'GridLineStyle',':');

