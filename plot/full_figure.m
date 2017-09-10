% Creat a figure occupying the full screen
function h_fig = full_figure(hold_flag)

scrsz = get(0,'ScreenSize');
h_fig = figure('Position',[10 40 scrsz(3)-20 scrsz(4)-110]);
if(~exist('hold_flag', 'var') || hold_flag)
    hold on;
end

