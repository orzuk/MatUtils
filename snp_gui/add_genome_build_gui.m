function varargout = add_genome_build_GUI(varargin)
% ADD_GENOME_BUILD_GUI M-file for add_genome_build_GUI.fig
%      ADD_GENOME_BUILD_GUI, by itself, creates a new ADD_GENOME_BUILD_GUI or raises the existing
%      singleton*.
%
%      H = ADD_GENOME_BUILD_GUI returns the handle to a new ADD_GENOME_BUILD_GUI or the handle to
%      the existing singleton*.
%
%      ADD_GENOME_BUILD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADD_GENOME_BUILD_GUI.M with the given input arguments.
%
%      ADD_GENOME_BUILD_GUI('Property','Value',...) creates a new ADD_GENOME_BUILD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before add_genome_build_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to add_genome_build_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help add_genome_build_GUI

% Last Modified by GUIDE v2.5 20-Sep-2007 14:09:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @add_genome_build_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @add_genome_build_GUI_OutputFcn, ...
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


% --- Executes just before add_genome_build_GUI is made visible.
function add_genome_build_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to add_genome_build_GUI (see VARARGIN)

% Choose default command line output for add_genome_build_GUI
handles.output = hObject;

load(fullfile('..','database','present_genome_builds.mat'));

builds_list=cell(size(genome_builds,1),1);
for i=1:size(genome_builds,1)
    builds_list{i}=[genome_builds{i,1} ' (' genome_builds{i,3} ',' genome_builds{i,2} ')'];   %hg17 (May 2004, NCBI Build 35)
end

set(handles.edit_supported_builds,'string',builds_list);

handles.current_genome_builds=genome_builds;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes add_genome_build_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = add_genome_build_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_supported_builds_Callback(hObject, eventdata, handles)
% hObject    handle to edit_supported_builds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_supported_builds as text
%        str2double(get(hObject,'String')) returns contents of edit_supported_builds as a double


% --- Executes during object creation, after setting all properties.
function edit_supported_builds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_supported_builds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_continue.
function pushbutton_continue_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_continue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

genome_assembly = handles.genome_build{1}; %'hg18';
new_DB_folder=fullfile('..','new_DB');
mkdir(new_DB_folder);

listbox_chips=get(handles.listbox_chips,'string');
selected_chips_idx=get(handles.listbox_chips,'value');

%1
create_chr_bands_data(genome_assembly, new_DB_folder)

%2
combine_refgene_reflink('refGene.txt','refLink_11cols.txt',genome_assembly);

ucsc_gene_loc_table_into_mat_with_exons(genome_assembly)


%3
% chip_name='Xba';
% annot_file='Mapping50K_Xba240.na22.annot.csv';

% chip_name='Hind';
% annot_file='Mapping50K_Hind240.na22.annot.csv';
% 
chip_name='Sty';
annot_file='Mapping250K_Sty.na22.annot.csv';
% 
% chip_name='Nsp';
% annot_file='Mapping250K_Nsp.na22.annot.csv';

load_chip_allele_info_table(chip_name, genome_assembly, annot_file);
build_snp_chip_annot_mat(genome_assembly, chip_name);

save(fullfile(new_DB_folder,'settings.mat'),'genome_assembly');

% --- Executes on selection change in listbox_chips.
function listbox_chips_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_chips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_chips contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_chips


% --- Executes during object creation, after setting all properties.
function listbox_chips_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_chips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_go_to_affy2.
function pushbutton_go_to_affy2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_go_to_affy2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web('http://www.affymetrix.com/support/index.affx','-browser')


function edit_add_build_Callback(hObject, eventdata, handles)
% hObject    handle to edit_add_build (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_add_build as text
%        str2double(get(hObject,'String')) returns contents of edit_add_build as a double


% --- Executes during object creation, after setting all properties.
function edit_add_build_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_add_build (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_add_build_ok.
function pushbutton_add_build_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_build_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

genome_build=get(handles.edit_add_build,'string');
idx = strfind(genome_build, ',');
genome_build_name = genome_build(1:idx-1);
genome_build_date = genome_build(idx+2:end);

if find(strcmp(handles.current_genome_builds(:,2), genome_build_name))
    button = questdlg([genome_build_name ' is already present in the database. Would you like to add it anyway?']);
    if ~strcmpi(button,'yes')
        return;
    end
end

web('http://genome.ucsc.edu/FAQ/FAQreleases#release1','-browser')

ucsc_version = inputdlg(['Please insert the UCSC Version for ' genome_build_name ...
    ' . This could be found in the just opened browser at http://genome.ucsc.edu/FAQ/FAQreleases#release1'],'UCSC Version');

if isempty(ucsc_version)
    return;
end
if isempty(ucsc_version{1})
    errordlg('Invalid UCSC version.','Error','modal');
    return;
end

genome_build=cell(1,3);
genome_build{1,1}= strtrim(ucsc_version{1});
genome_build{1,2}= genome_build_name;
genome_build{1,3}= genome_build_date;

handles.genome_build = genome_build;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_go_to_affy1.
function pushbutton_go_to_affy1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_go_to_affy1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web('http://www.affymetrix.com/support/index.affx','-browser')


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

