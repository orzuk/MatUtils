function varargout = display_raw_data_GUI(varargin)
% DISPLAY_RAW_DATA_GUI M-file for display_raw_data_GUI.fig
%      DISPLAY_RAW_DATA_GUI, by itself, creates a new DISPLAY_RAW_DATA_GUI or raises the existing
%      singleton*.
%
%      H = DISPLAY_RAW_DATA_GUI returns the handle to a new DISPLAY_RAW_DATA_GUI or the handle to
%      the existing singleton*.
%
%      DISPLAY_RAW_DATA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAY_RAW_DATA_GUI.M with the given input arguments.
%
%      DISPLAY_RAW_DATA_GUI('Property','Value',...) creates a new DISPLAY_RAW_DATA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before display_raw_data_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to display_raw_data_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help display_raw_data_GUI

% Last Modified by GUIDE v2.5 22-Aug-2007 12:12:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @display_raw_data_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @display_raw_data_GUI_OutputFcn, ...
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


% --- Executes just before display_raw_data_GUI is made visible.
function display_raw_data_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to display_raw_data_GUI (see VARARGIN)

% Choose default command line output for display_raw_data_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes display_raw_data_GUI wait for user response (see UIRESUME)
uiwait(handles.display_raw_data_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = display_raw_data_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_start_Callback(hObject, eventdata, handles)
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_start as text
%        str2double(get(hObject,'String')) returns contents of edit_start as a double


% --- Executes during object creation, after setting all properties.
function edit_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_end_Callback(hObject, eventdata, handles)
% hObject    handle to edit_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_end as text
%        str2double(get(hObject,'String')) returns contents of edit_end as a double


% --- Executes during object creation, after setting all properties.
function edit_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_current_region.
function radiobutton_current_region_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_current_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_current_region


% --- Executes on button press in radiobutton_choose_region.
function radiobutton_choose_region_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_choose_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_choose_region


% --- Executes on button press in checkbox_disease.
function checkbox_disease_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_disease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_disease


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radiobutton_current_region,'value');
    params.region='current';
else
    params.region='choose';
    params.chrom=get(handles.popupmenu_chr,'value');
    params.start=str2num(get(handles.edit_start,'string'));
    params.end=str2num(get(handles.edit_end,'string'));
    if isempty(params.start) || params.start<0 || isempty(params.end) || params.end<0 || params.start >= params.end
        errordlg('Region Start and End must be non-negative. End must be greater than Start.','Invalid Region','modal');
        return;
    end
end
    

params.disease=get(handles.checkbox_disease,'value');
params.normal=get(handles.checkbox_normal,'value');
if ~(params.disease || params.normal)
    errordlg('At least Show Disease or Normal Samples must be checked.','Invalid Selection','modal');
    return;
end



params.hmm_data=get(handles.checkbox_hmm_data,'value');

out1.params=params;
out1.figure1=handles.display_raw_data_GUI;
out1.ok=1;


handles.output=out1;
guidata(hObject,handles);
uiresume;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out1.figure1=handles.display_raw_data_GUI;
out1.ok=0;
handles.output=out1;
guidata(hObject,handles);
uiresume;


% --- Executes on selection change in popupmenu_chr.
function popupmenu_chr_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_chr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chr


% --- Executes during object creation, after setting all properties.
function popupmenu_chr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close display_raw_data_GUI.
function display_raw_data_GUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to display_raw_data_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);

out1.figure1=handles.display_raw_data_GUI;
out1.ok=0;
handles.output=out1;
guidata(hObject,handles);
uiresume;




% --- Executes on button press in checkbox_hmm_data.
function checkbox_hmm_data_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_hmm_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_hmm_data




% --- Executes on button press in checkbox_normal.
function checkbox_normal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_normal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_normal


