function varargout = com_aber_options(varargin)
% COM_ABER_OPTIONS M-file for com_aber_options.fig
%      COM_ABER_OPTIONS, by itself, creates a new COM_ABER_OPTIONS or raises the existing
%      singleton*.
%
%      H = COM_ABER_OPTIONS returns the handle to a new COM_ABER_OPTIONS or the handle to
%      the existing singleton*.
%
%      COM_ABER_OPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COM_ABER_OPTIONS.M with the given input arguments.
%
%      COM_ABER_OPTIONS('Property','Value',...) creates a new COM_ABER_OPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before com_aber_options_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to com_aber_options_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help com_aber_options

% Last Modified by GUIDE v2.5 07-Aug-2007 16:12:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @com_aber_options_OpeningFcn, ...
                   'gui_OutputFcn',  @com_aber_options_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before com_aber_options is made visible.
function com_aber_options_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to com_aber_options (see VARARGIN)

% Choose default command line output for com_aber_options
handles.output = hObject;

params=varargin{1};

% if params.current.del==1
set(handles.checkbox_del,'value',params.current.del);
% else
%     set(handles.checkbox_del,'value',0);
% end

% if params.current.run_chips_together==1
set(handles.checkbox_run_chips_together,'value',params.current.run_chips_together);
% else
%     set(handles.checkbox_run_chips_together,'value',0);
% end

set(handles.edit_del_thresh,'string',num2str(params.current.del_thresh));
set(handles.edit_amp_thresh,'string',num2str(params.current.amp_thresh));
set(handles.edit_fdr,'string',num2str(params.current.fdr));
set(handles.edit_min_stretch,'string',num2str(params.current.min_stretch));
set(handles.edit_max_stretch,'string',num2str(params.current.max_stretch));

 
set(handles.checkbox_disease,'value',params.current.disease);
set(handles.checkbox_normals,'value',params.current.normals);
set(handles.checkbox_remove_del_arm,'value',params.current.remove_del_arm);
set(handles.checkbox_remove_amp_arm,'value',params.current.remove_amp_arm);
set(handles.edit_fraction_arm,'string',num2str(params.current.fraction_arm));
set(handles.edit_webpage_ucsc_input,'string',params.current.ucsc_webpage);
set(handles.edit_mounted_folder,'string',params.current.ucsc_mounted_folder);

handles.params_default=params.default;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes com_aber_options wait for user response (see UIRESUME)
uiwait(handles.com_aber_options);


% --- Outputs from this function are returned to the command line.
function varargout = com_aber_options_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_del.
function checkbox_del_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_del


% --- Executes on button press in checkbox_run_chips_together.
function checkbox_run_chips_together_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_chips_together (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_chips_together



function edit_fdr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fdr as text
%        str2double(get(hObject,'String')) returns contents of edit_fdr as a double


% --- Executes during object creation, after setting all properties.
function edit_fdr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

params.del=get(handles.checkbox_del,'value');
params.run_chips_together=get(handles.checkbox_run_chips_together,'value');

% if params.del==0 && params.amp==0
%     errordlg('At least deletion or amplification must be selected','Invalid Selection','modal');
%     return;
% end

params.del_thresh=str2num(get(handles.edit_del_thresh,'string'));
params.amp_thresh=str2num(get(handles.edit_amp_thresh,'string'));

if isempty(params.del_thresh) || params.del_thresh<0 || isempty(params.amp_thresh) || params.amp_thresh<0
    errordlg('Deletion and amplification thresholds must be non-negative numbers','Invalid Thresholds','modal');
    return;
end

params.fdr=str2num(get(handles.edit_fdr,'string'));
if isempty(params.fdr) || params.fdr<0 || params.fdr>1
    errordlg('FDR Q must be a number between 0 to 1','Invalid FDR Q','modal');
    return;
end

params.min_stretch=str2num(get(handles.edit_min_stretch,'string'));
params.max_stretch=str2num(get(handles.edit_max_stretch,'string'));
if isempty(params.min_stretch) || params.min_stretch<=0 || mod(params.min_stretch,1)~=0 || ...
    isempty(params.max_stretch) || params.max_stretch<=0 || mod(params.max_stretch,1)~=0 || params.min_stretch>params.max_stretch
    errordlg('Min and Max stretches must be a valid range of integers','Invalid Stretch','modal');
    return;
end

 
params.disease=get(handles.checkbox_disease,'value');
params.normals=get(handles.checkbox_normals,'value');

if ~params.disease && ~params.normals
    errordlg('At least "run on disease" or "run on normals" must be checked','Error','modal');
    return;
end

params.remove_del_arm=get(handles.checkbox_remove_del_arm,'value');
params.remove_amp_arm=get(handles.checkbox_remove_amp_arm,'value');

params.fraction_arm=str2num(get(handles.edit_fraction_arm,'string'));
if isempty(params.fraction_arm) || params.fraction_arm<0 || params.fraction_arm>1
    errordlg('Fraction Whole Arm must be a number between 0 to 1','Invalid Fraction Whole Arm','modal');
    return;
end

params.ucsc_webpage=strtrim(get(handles.edit_webpage_ucsc_input,'string'));
params.ucsc_mounted_folder=strtrim(get(handles.edit_mounted_folder,'string'));
if ~isempty(params.ucsc_mounted_folder) && ~exist(params.ucsc_mounted_folder,'dir')
    errordlg('Pleae enter a valid mounted folder or leave empty for using the output folder.','Invalid Mounted Folder','modal');
    return;
end

out1.params=params;
out1.figure1=handles.com_aber_options;
out1.ok=1;


handles.output=out1;
guidata(hObject,handles);
uiresume;




% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out1.figure1=handles.com_aber_options;
out1.ok=0;
handles.output=out1;
guidata(hObject,handles);
uiresume;





% --- Executes when user attempts to close com_aber_options.
function com_aber_options_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to com_aber_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);

out1.figure1=handles.com_aber_options;
out1.ok=0;
handles.output=out1;
guidata(hObject,handles);
uiresume;




% --- Executes on button press in pushbutton_defaults.
function pushbutton_defaults_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.params_default.del==1
    set(handles.checkbox_del,'value',1);
else
    set(handles.checkbox_del,'value',0);
end

if handles.params_default.run_chips_together==1
    set(handles.checkbox_run_chips_together,'value',1);
else
    set(handles.checkbox_run_chips_together,'value',0);
end


set(handles.edit_del_thresh,'string',num2str(handles.params_default.del_thresh));
set(handles.edit_amp_thresh,'string',num2str(handles.params_default.amp_thresh));
set(handles.edit_fdr,'string',num2str(handles.params_default.fdr));
set(handles.edit_min_stretch,'string',num2str(handles.params_default.min_stretch));
set(handles.edit_max_stretch,'string',num2str(handles.params_default.max_stretch));

set(handles.checkbox_disease,'value',handles.params_default.disease);
set(handles.checkbox_normals,'value',handles.params_default.normals);

set(handles.checkbox_remove_del_arm,'value',handles.params_default.remove_del_arm);
set(handles.checkbox_remove_amp_arm,'value',handles.params_default.remove_amp_arm);
set(handles.edit_fraction_arm,'string',num2str(handles.params_default.fraction_arm));
set(handles.edit_webpage_ucsc_input,'string',handles.params_default.ucsc_webpage);
set(handles.edit_mounted_folder,'string',handles.params_default.ucsc_mounted_folder);



function edit_amp_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_amp_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_amp_thresh as text
%        str2double(get(hObject,'String')) returns contents of edit_amp_thresh as a double


% --- Executes during object creation, after setting all properties.
function edit_amp_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_amp_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_del_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_del_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_del_thresh as text
%        str2double(get(hObject,'String')) returns contents of edit_del_thresh as a double


% --- Executes during object creation, after setting all properties.
function edit_del_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_del_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_webpage_ucsc_input_Callback(hObject, eventdata, handles)
% hObject    handle to edit_webpage_ucsc_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_webpage_ucsc_input as text
%        str2double(get(hObject,'String')) returns contents of edit_webpage_ucsc_input as a double


% --- Executes during object creation, after setting all properties.
function edit_webpage_ucsc_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_webpage_ucsc_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mounted_folder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mounted_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mounted_folder as text
%        str2double(get(hObject,'String')) returns contents of edit_mounted_folder as a double


% --- Executes during object creation, after setting all properties.
function edit_mounted_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mounted_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_max_stretch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_stretch as text
%        str2double(get(hObject,'String')) returns contents of edit_max_stretch as a double


% --- Executes during object creation, after setting all properties.
function edit_max_stretch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_stretch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_stretch as text
%        str2double(get(hObject,'String')) returns contents of edit_min_stretch as a double


% --- Executes during object creation, after setting all properties.
function edit_min_stretch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_remove_del_arm.
function checkbox_remove_del_arm_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_remove_del_arm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_remove_del_arm


% --- Executes on button press in checkbox_remove_amp_arm.
function checkbox_remove_amp_arm_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_remove_amp_arm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_remove_amp_arm



function edit_fraction_arm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fraction_arm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fraction_arm as text
%        str2double(get(hObject,'String')) returns contents of edit_fraction_arm as a double


% --- Executes during object creation, after setting all properties.
function edit_fraction_arm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fraction_arm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_remove_del_arm.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_remove_del_arm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_remove_del_arm


% --- Executes on button press in checkbox_remove_amp_arm.
function checkbox_remove_amp_arms_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_remove_amp_arm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_remove_amp_arm


% --- Executes on button press in checkbox_del.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_del




% --- Executes on button press in checkbox_disease.
function checkbox_disease_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_disease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_disease


% --- Executes on button press in checkbox_normals.
function checkbox_normals_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_normals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_normals


