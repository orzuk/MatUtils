% --- Executes on mouse button press on axes
function mapAxesClick_Callback(hObject, eventdata)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(gcbo);

% remove old rect
f = findobj(handles.subplots(get(handles.chromAxes, 'UserData')),'Tag','choice_rect');
if ~isempty(f)
    delete(f);
end

% set(hObject, 'Selected', 'on', 'SelectionHighlight', 'on');
chrom = get(gca, 'UserData');
cla(handles.chromAxes);
set(handles.uipanel_choose_view, 'SelectedObject', handles.radioBands);
handles.h_lines_standAlone=drawChrom(handles, chrom, 1, 0);
% h=get(handles.chromAxes, 'Children');
set(handles.chromAxes, 'ButtonDownFcn', 'chromAxesClick_Callback');
guidata(gcbo, handles);

