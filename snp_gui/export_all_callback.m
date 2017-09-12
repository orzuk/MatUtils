function export_all_Callback(hObject, eventdata)
handles=guidata(gcbo);

figure      
% drawChromMap(handles, 0);
% set(get(gcf, 'children'), 'ButtonDownFcn', []);


left=0.05;
bottom=0.25;
width=0.9;
height=0.73;

h = [];

subplots=handles.subplots;

for i=1:length(subplots)
    h(i) = axes('Position',[left+0.02 bottom+0.01+(1/(length(subplots)+1))*height*(i-1) width-left-0.02 (1/(length(subplots)+1))*height*0.9]);
    set(gca, 'Color', 'none');
    set(gca,'XTickLabel',[], 'YTickLabel',[]);
    set(gca,'XTick',[], 'YTick',[]);
    set(gca, 'position', get(subplots(i), 'position'), 'xLim', get(subplots(i), 'xLim'), 'yLim', get(subplots(i), 'yLim'));
    axis off; 
    h=get(subplots(i), 'Children');
    new_handle = copyobj(h,gca);
    set(new_handle, 'ButtonDownFcn', []);
end
