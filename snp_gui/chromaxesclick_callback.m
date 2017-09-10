% --- Executes on button press in chromAxes
function chromAxesClick_Callback(hObject, eventdata)
% hObject    handle to zoomOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbf);
point1 = get(handles.chromAxes,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(handles.chromAxes,'CurrentPoint');    % button up detected
chrom = get(handles.chromAxes, 'UserData');

point = [point1(1,1:2), point2(1,1:2)];
xmin = min(point(1), point(3));
xmax = max(point(1), point(3));
zoomIntoChrom(handles,xmin,xmax,chrom);
guidata(gcbo, handles);



