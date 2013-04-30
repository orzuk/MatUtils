function zoomIntoChrom(handles,xmin,xmax,chrom)


if xmax-xmin<1 %==0
    return
end

% remove old rect
f = findobj(handles.subplots(chrom),'Tag','choice_rect');
if ~isempty(f)
    delete(f);
end


axes(handles.subplots(chrom));
ax=axis;
y = ax(3);
height = ax(4)-ax(3);

xmin = max(xmin, ax(1));
xmax = min(xmax, max(handles.chr_band_end(handles.chr_num==chrom)));
width = xmax-xmin;

p = rectangle('Position',[xmin y width height], 'LineWidth', 2, 'EdgeColor', 'g', 'Tag', 'choice_rect');
 
axes(handles.chromAxes);

%update slider steps
set(handles.chromAxes, 'XLim', [xmin xmax]);
set(handles.chromSlider, 'Value', xmin);
bigStep = min((xmax-xmin)/(get(handles.chromSlider,'Max') - get(handles.chromSlider,'Min')), 0.1);
smallStep = min((xmax-xmin)/(5*(get(handles.chromSlider,'Max') - get(handles.chromSlider,'Min'))), 00.1);
set(handles.chromSlider, 'SliderStep', [smallStep, bigStep]);

updateChromView(handles); %was added on 02/01/2007

