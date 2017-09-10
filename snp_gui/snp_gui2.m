function varargout = SNP_GUI(varargin)
% SNP_GUI M-file for SNP_GUI.fig
%      SNP_GUI, by itself, creates a new SNP_GUI or raises the existing
%      singleton*.
%
%      H = SNP_GUI returns the handle to a new SNP_GUI or the handle to
%      the existing singleton*.
%
%      SNP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNP_GUI.M with the given input arguments.
%
%      SNP_GUI('Property','Value',...) creates a new SNP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SNP_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SNP_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's tools_item_find_gene menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SNP_GUI

% Last Modified by GUIDE v2.5 27-Jun-2007 17:30:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SNP_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SNP_GUI_OutputFcn, ...
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


% --- Executes just before SNP_GUI is made visible.
function SNP_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SNP_GUI (see VARARGIN)

% Choose default command line output for SNP_GUI
handles.output = hObject;

%handles.data = [];
%handles.chr = [];
%handles.location = [];
% handles.start_over = [];
% handles.end_over = [];
% handles.start_under = [];
% handles.end_under = [];
% handles.qvals_over = [];
% handles.qvals_under = [];
handles.subplots = [];

% menu = uimenu('Label', 'Edit');
% item1 = uimenu(menu, 'Label', 'Export', 'Callback', 'export_all_Callback');

genome_assembly = get_genome_assembly();
load(fullfile('..','database',['chr_data_' genome_assembly '.mat']));
handles.chr_band_end = chr_band_end;
handles.chr_band_start = chr_band_start;
handles.chr_band_names = chr_band_names;
handles.chr_num = chr_num;
handles.end_p_location = end_p_location;
load(fullfile('..','database',['refgenes_' genome_assembly '_with_exons.mat']));
handles.annot.symbols = gene_symbols;
handles.annot.chr = chr;
handles.annot.start = loc_start;
handles.annot.end = loc_end;
% load(fullfile('data','snp_symbols.mat'));
% char_symbols=char(symbols);
% symbols=cellstr(char_symbols(:,2:end));
% handles.symbols = symbols;

set(handles.uipanel_choose_view, 'UserData', 'location');
set(handles.uipanel_choose_view, 'SelectedObject', handles.radioBands);

[handles.subplots handles.h_lines]=drawChromMap(handles, 1);

%default paramters for find common aberrations
handles.com_aber_params_default.del=1;
handles.com_aber_params_default.run_chips_together=0;
handles.com_aber_params_default.del_thresh=1.6;
handles.com_aber_params_default.amp_thresh=2.5;
handles.com_aber_params_default.fdr=0.05;
handles.com_aber_params_default.ucsc_webpage='http://www.weizmann.ac.il/';
handles.com_aber_params_default.ucsc_mounted_folder='';

handles.com_aber_params=handles.com_aber_params_default;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SNP_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SNP_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_load_norm.
function pushbutton_load_norm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chip_list=get(handles.popupmenuChipName,'string');
chip_name=chip_list{get(handles.popupmenuChipName,'value')};
if strcmp(chip_name,'Hind50k ( from Affy 100K)')
    chip_name='Hind';
elseif strcmp(chip_name,'Xba50k ( from Affy 100K)')
    chip_name='Xba';
else
    errordlg('Invalid chip','Chip Error','modal');
    return;
end

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

samplesFile=deblank(get(handles.edit_samp_file,'string'));
if ~exist(samplesFile,'file')
    errordlg('Samples file not found','File Error','modal');
    return;
end

mat=loadCellFile(samplesFile);
if size(mat,2)<3 || ~strcmpi(mat{1,1},'sample') || ~strcmpi(mat{1,2},'cel file')  || ~strcmpi(mat{1,3},'gender')
    errordlg({'First column title must be "sample"' 'Second column title must be "CEL file"' 'Third column title must be "Gender"'},'File Error','modal');
    return;
end
    
samples=mat(2:end,1);
arrays=mat(2:end,2);
genders=mat(2:end,3);

for i=1:length(arrays)
    if isempty(fileparts(arrays{i}))
        arrays{i}=fullfile(cd,'..','data',arrays{i});
    end
    if ~exist(arrays{i},'file')
        errordlg(['CEL file: ' arrays{i} ' not found'],'File Error','modal');
        return;
    end
end

button = questdlg('Would you like to run HMM right after loading?','run HMM','Yes','No','Yes');
if isempty(button)
    return;
end
if strcmp(button,'Yes')
    answer=inputdlg({'Number of EM iterations' 'Number of EM starting points' 'EM_tolerance' 'Use genotype correlations (y or n)'} ...
        ,'HMM Par.',1,...
        {'25' '5' '0.000001' 'n'});

    if isempty(answer)
        return;
    end
    handles.HMMParamsStruct.num_EM_iters = str2num(answer{1}); %50;
    handles.HMMParamsStruct.num_EM_starting_points = str2num(answer{2}); %10;
    handles.HMMParamsStruct.EM_tolerance = str2num(answer{3}); %0.000001;
    if strcmpi(answer{4},'y')
        handles.HMMParamsStruct.UseGenotypesCorrs=1;
    elseif strcmpi(answer{4},'n')
        handles.HMMParamsStruct.UseGenotypesCorrs=0;
    else
        errordlg('"Use genotype correlations" field must be either y or n','Parameter Error','modal');
        return;
    end
    handles.run_HMM_after_load=1;
else
    handles.run_HMM_after_load=0;
end


normal_ind = get_normal_copy_num_samples_ind(samples); female_ind = strmatch('F', genders); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_cel_files_flag=1; prep_block_flag=1; % Can change to zero after first time ..
RLMM_Normalize(workPath, samples, normal_ind,  female_ind, arrays, chip_name, load_cel_files_flag, prep_block_flag); % Do not give correct dir (it will overide files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diploid_ind = get_diploid_ind(samples);
%load_and_norm_data(samples, arrays, workPath, diploid_ind, chip_name);

if handles.run_HMM_after_load
    disp('Loading complete. Running HMM ...');
    pushbutton_run_hmm_Callback(handles.pushbutton_run_hmm, eventdata, handles)
else
    msgbox('Load and Normalization have finished.','Load and Normalization Finished','modal');
end


% --- Executes on selection change in popupmenuChipName.
function popupmenuChipName_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuChipName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuChipName contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuChipName


% --- Executes during object creation, after setting all properties.
function popupmenuChipName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuChipName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_run_hmm.
function pushbutton_run_hmm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run_hmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

genome_assembly = get_genome_assembly();
chip_list=get(handles.popupmenuChipName,'string');
chip_name=chip_list{get(handles.popupmenuChipName,'value')};

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

samplesFile=deblank(get(handles.edit_samp_file,'string'));
if exist(samplesFile)~=2
    errordlg('Samples file not found','File Error','modal');
    return;
end

samples=loadCellFile(samplesFile);
if ~strcmpi(samples{1,1},'sample')
    errordlg('First column title must be "sample"','File Error','modal');
    return;
end

if isfield(handles,'run_HMM_after_load') && handles.run_HMM_after_load
    HMMParamsStruct=handles.HMMParamsStruct;
else
    answer=inputdlg({'Number of EM iterations' 'Number of EM starting points' 'EM_tolerance' 'Use genotype correlations (y or n)'} ...
        ,'HMM Par.',1,...
        {'25' '5' '0.000001' 'n'});

    if isempty(answer)
        return;
    end
    HMMParamsStruct.num_EM_iters = str2num(answer{1}); %50;
    HMMParamsStruct.num_EM_starting_points = str2num(answer{2}); %10;
    HMMParamsStruct.EM_tolerance = str2num(answer{3}); %0.000001;
    if strcmpi(answer{4},'y')
        HMMParamsStruct.UseGenotypesCorrs=1;
    elseif strcmpi(answer{4},'n')
        HMMParamsStruct.UseGenotypesCorrs=0;
    else
        errordlg('"Use genotype correlations" field must be either y or n','Parameter Error','modal');
        return;
    end
end

samples=samples(2:end,1);

if strcmp(chip_name,'Hind50k ( from Affy 100K)')
    LDStruct=load(fullfile('..','database','LD_Hind.mat'));
    handles.chip_annot=load(fullfile('..','database',['Hind_annot_data_' genome_assembly '.mat']));
    handles.ChipName='Hind';
elseif strcmp(chip_name,'Xba50k ( from Affy 100K)')
    LDStruct=load(fullfile('..','database','LD_xba.mat'));
    handles.chip_annot=load(fullfile('..','database',['xba_annot_data_' genome_assembly '.mat']));
    handles.ChipName='Xba';
else
    errordlg('Invalid chip','Chip Error','modal');
    return;
end
    
HMMParamsStruct.x_dim = 3; HMMParamsStruct.y_dim = 1;


HMMParamsStruct.learn_model_EM = 1;
HMMParamsStruct.ChromosomesToRun = 1:22 ;
HMMParamsStruct.ChromosomesToPlot = 1:22 ; % Plot

HMMParamsStruct.derich = 1.0/120.0;

    
%% OutputFilesNames = FakeAfterRunHMMFile(workPath, samples, LDStruct, handles.chip_annot, HMMParamsStruct);
OutputFilesNames = LearnHMMChromFromData(workPath, samples, LDStruct, handles.chip_annot, HMMParamsStruct);

% handles.OutputFilesNames=OutputFilesNames;
guidata(hObject, handles);

save(fullfile(workPath,['after_run_hmm_' handles.ChipName '.mat']),'OutputFilesNames','samples');

set(handles.popupmenu_samples,'string',['Average'; samples]);

msgbox('Run HMM has finished.','Run HMM Finished','modal');


% --- Executes on button press in pushbutton_view_res.
function pushbutton_view_res_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_view_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

genome_assembly = get_genome_assembly();
workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

if get(handles.checkbox_all_samples','value')
    if ~exist(fullfile(workPath,'output','stretches_disp.mat'),'file')
        errordlg('Find Common Aberrations not been performed','Results Error','modal');
        return;
    end
    load(fullfile(workPath,'output','stretches_disp.mat'));
else
    chip_list=get(handles.popupmenuChipName,'string');
    chip_name=chip_list{get(handles.popupmenuChipName,'value')};
    if strcmp(chip_name,'Hind50k ( from Affy 100K)')
        chip_name='Hind';
    elseif strcmp(chip_name,'Xba50k ( from Affy 100K)')
        chip_name='Xba';
    else
        errordlg('Invalid chip','Chip Error','modal');
        return;
    end

    file_name=fullfile(workPath,['after_run_hmm_' chip_name '.mat']);
    if ~exist(file_name,'file')
        errordlg(['Run HMM has not been performed. ' file_name ' not found.'],'Results Error','modal');
        return;
    end
    load(file_name);
    if length(OutputFilesNames.disp_files)~=length(samples)+1
        errordlg('length(OutputFilesNames.disp_files) must be length(samples)+1','Samples Error','modal');
        return;
    end

    samples_from_list=get(handles.popupmenu_samples,'string');
    samples_from_list(1)=[];
    if ~isequal(samples_from_list,samples)
        set(handles.popupmenu_samples,'string',['Average'; samples]);
        msgbox('Please choose a sample and click "View Results" again.','View Results','modal');
        return;
    end
    
    %load(OutputFilesNames.disp_files{get(handles.popupmenu_samples,'value')});
    load(fullfile(workPath,OutputFilesNames.disp_files{get(handles.popupmenu_samples,'value')}));
end

if strcmpi(DispStruct.ChipName,'Hind')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Hind')
        handles.chip_annot=load(fullfile('..','database',['Hind_annot_data_' genome_assembly '.mat']));
        handles.ChipName='Hind';
    end
elseif strcmpi(DispStruct.ChipName,'Xba')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Xba')
        handles.chip_annot=load(fullfile('..','database',['Xba_annot_data_' genome_assembly '.mat']));
        handles.ChipName='Xba';
    end
elseif strcmp(DispStruct.ChipName,'Hind_and_Xba')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Hind_and_Xba')
        handles.chip_annot=load_Hind_and_Xba;
        handles.ChipName='Hind_and_Xba';
    end
else
    errordlg(['Invalid chip ' DispStruct.ChipName],'Chip Error','modal');
    return;
end

%load data\SampleChromsDataStructForDisplay.mat;
handles.GenomeBuild=DispStruct.GenomeBuild;
handles.SampleName=DispStruct.SampleName;
handles.ChromsData=DispStruct.Chrom;
% handles.data = data;
% handles.chr = chr;
% handles.location = location;
% handles.start_over = start_over;
% handles.end_over = end_over;
% handles.qvals_over = qvals_over;
% handles.start_under = start_under;
% handles.end_under = end_under;
% handles.qvals_under = qvals_under;
%assume that data is after log2 and smoothing
[handles.subplots handles.h_lines handles.min_intensity_seg_all handles.max_intensity_seg_all]=drawChromMap(handles, 0);
%handles.subplots = h;
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_samples.
function popupmenu_samples_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_samples contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_samples


% --- Executes during object creation, after setting all properties.
function popupmenu_samples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on slider movement.
function chromSlider_Callback(hObject, eventdata, handles)
% hObject    handle to chromSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

axes(handles.chromAxes);

chrom = get(handles.chromAxes, 'UserData');
% remove old rect
f = findobj(handles.subplots(chrom),'Tag','choice_rect');
if ~isempty(f)
    delete(f);
end

steps = get(hObject,'SliderStep');
limits = get(gca,'Xlim');
width = steps(2)*(get(hObject,'Max') - get(hObject,'Min'));
if (get(hObject,'Max')==limits(2) & get(hObject,'Min')==limits(1)) %whole chromosome view
    set(hObject,'Value', limits(1));
    return;
end


width = limits(2)- limits(1);
if get(hObject,'Value')<limits(1)
    start_pos = get(hObject,'Value');
    end_pos = get(hObject,'Value') + width;
elseif get(hObject,'Value')>limits(2)
    start_pos = get(hObject,'Value') - width;
    end_pos = get(hObject,'Value') ;
else
    return;
end

set(gca, 'xLim', [start_pos end_pos]);

axes(handles.subplots(chrom));
ax=axis;
y = ax(3);
height = ax(4)-ax(3);
p = rectangle('Position',[start_pos y width height], 'LineWidth', 2, 'EdgeColor', 'g','Tag','choice_rect');

axes(handles.chromAxes);
guidata(hObject, handles);
updateChromView(handles);


% --- Executes during object creation, after setting all properties.
function chromSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chromSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --------------------------------------------------------------------
function uipanel_choose_view_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_choose_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateChromView(handles);




% --------------------------------------------------------------------
function menu_export_Callback(hObject, eventdata, handles)
% hObject    handle to menu_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export_item_export_all_Callback(hObject, eventdata, handles)
% hObject    handle to export_item_export_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

export_all_Callback;


% --------------------------------------------------------------------
function export_item_export_region_Callback(hObject, eventdata, handles)
% hObject    handle to export_item_export_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


figure
copyobj(handles.chromAxes,gcf);


% --------------------------------------------------------------------
function tools_item_find_gene_Callback(hObject, eventdata, handles)
% hObject    handle to tools_item_find_gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'chip_annot')
    errordlg('Please load data','Data not Loaded','modal');
    return;
end
    

handles = guidata(hObject);
geneSymbol = deblank(inputdlg('Enter Gene Symbol'));
if isempty(geneSymbol)
   return;
end

snpIx=strmatch(lower(geneSymbol), lower(handles.chip_annot.snps_gene_symbols), 'exact');
geneIx=strmatch(lower(geneSymbol), lower(handles.annot.symbols), 'exact');
if isempty(snpIx) || isempty(geneIx)
    errordlg(['Cannot find gene symbol: ' geneSymbol],'Gene not Found','modal');
    return;
end
xmin=min(handles.chip_annot.chr_loc_vec(min(snpIx))-100000, handles.annot.start(geneIx));
xmax=max(handles.chip_annot.chr_loc_vec(max(snpIx))+100000, handles.annot.end(geneIx));
chrom=handles.chip_annot.chr_vec(snpIx(1));

handles.h_lines_standAlone=drawChrom(handles, chrom, 1, 0);
zoomIntoChrom(handles,xmin,xmax,chrom);
guidata(gcbo, handles);


% --------------------------------------------------------------------
function tools_item_genecards_Callback(hObject, eventdata, handles)
% hObject    handle to tools_item_genecards (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

geneSymbol = deblank(inputdlg('Enter Gene Symbol'));
if isempty(geneSymbol)
   return;
end

web(['http://www.genecards.org/cgi-bin/carddisp.pl?gene=' geneSymbol{1}],'-browser')

% --------------------------------------------------------------------
function tools_item_ucsc_Callback(hObject, eventdata, handles)
% hObject    handle to tools_item_ucsc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'GenomeBuild')
    errordlg('Data should be loaded.','Genome Build Error','modal');
    return;
end

chrom=num2str(get(handles.chromAxes,'UserData'));
if isempty(chrom)
    errordlg('A chromosome should be chosen.','Chromosome Error','modal');
    return;
end
    
xlim1=get(handles.chromAxes,'xlim');


web(['http://genome.ucsc.edu/cgi-bin/hgTracks?db=' lower(handles.GenomeBuild) ...
    '&position=chr' chrom ':' num2str(floor(xlim1(1))) '-' num2str(ceil(xlim1(2)))],'-browser');

% --------------------------------------------------------------------
function tools_menu_Callback(hObject, eventdata, handles)
% hObject    handle to tools_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes on button press in pushbuttonWorkPath.
function pushbuttonWorkPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonWorkPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------
cd('..'); % because of a bug in uigetdir 5/12/2006
directory_name = uigetdir(fullfile(cd,'data'),'Select Working Directory');
cd('src');
%-------------------------
if directory_name~=0
    set(handles.editWorkPath,'string',directory_name);
end

function editWorkPath_Callback(hObject, eventdata, handles)
% hObject    handle to editWorkPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWorkPath as text
%        str2double(get(hObject,'String')) returns contents of editWorkPath as a double


% --- Executes during object creation, after setting all properties.
function editWorkPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWorkPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_samp_file.
function pushbutton_samp_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_samp_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

files='*.txt';
workDir=get(handles.editWorkPath,'string');
if exist(workDir)==7
    files=fullfile(workDir,files);
end
[FileName,PathName] = uigetfile(files, 'Select Samples File');
if FileName~=0
    set(handles.edit_samp_file,'string',fullfile(PathName,FileName));
end




% --------------------------------------------------------------------
function options_menu_Callback(hObject, eventdata, handles)
% hObject    handle to options_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function com_aber_item_Callback(hObject, eventdata, handles)
% hObject    handle to com_aber_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

com_aber_params.current=handles.com_aber_params;
com_aber_params.default=handles.com_aber_params_default;

options=com_aber_options(com_aber_params);
delete(options.figure1);

if options.ok==1
    handles.com_aber_params=options.params;
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_find_com_aber.
function pushbutton_find_com_aber_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_find_com_aber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chip_list=get(handles.popupmenuChipName,'string');
chip_name=chip_list{get(handles.popupmenuChipName,'value')};
genome_assembly = get_genome_assembly();

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

samplesFile=deblank(get(handles.edit_samp_file,'string'));
if exist(samplesFile)~=2
    errordlg('Samples file not found','File Error','modal');
    return;
end

samples=loadCellFile(samplesFile);
if ~strcmpi(samples{1,1},'sample')
    errordlg('First column title must be "sample"','File Error','modal');
    return;
end
    
samples=samples(2:end,1);

if handles.com_aber_params.run_chips_together
    chip_name='Hind_and_Xba';
end

if strcmp(chip_name,'Hind50k ( from Affy 100K)')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Hind')
        handles.chip_annot=load(fullfile('..','database',['Hind_annot_data_' genome_assembly '.mat']));
        handles.ChipName='Hind';
    end
elseif strcmp(chip_name,'Xba50k ( from Affy 100K)')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Xba')
        handles.chip_annot=load(fullfile('..','database',['Xba_annot_data_' genome_assembly '.mat']));
        handles.ChipName='Xba';
    end
elseif strcmp(chip_name,'Hind_and_Xba')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Hind_and_Xba')
        handles.chip_annot=load_Hind_and_Xba;
        handles.ChipName='Hind_and_Xba';
    end
else
    errordlg('Invalid chip','Chip Error','modal');
    return;
end

guidata(hObject, handles);

AssignAllGlobalConstants();  

genes_db_f = ['refgenes_' genome_assembly '_with_exons.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load SNP annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snp_ids = handles.chip_annot.snp_ids;
chr_loc_vec = handles.chip_annot.chr_loc_vec;
chr_vec = handles.chip_annot.chr_vec;
snp_gene_symbols = handles.chip_annot.snps_gene_symbols;
%snps_gene_dist = handles.chip_annot.snps_gene_dist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load genes description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genes_db_struct = load(['../Database/' genes_db_f]); 
gene_descr = genes_db_struct.gene_descr;
gene_symbols = genes_db_struct.gene_symbols;
% make gene_descr for snp_gene_symbols
snp_gene_descr = snp_gene_symbols_into_descr(snp_gene_symbols, gene_symbols, gene_descr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run function that find stretches and saves output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%samples = {'TEL74_5_n';'TEL74_5_d';'TEL76_7_n';'TEL76_7_d'};

%     ;'HD78_9_n';'HD78_9_d';...
%     'TEL80_1_n';'TEL80_1_d';'HD82_3_n';'HD82_3_d';'TEL84_5_n';'TEL84_5_d';'HD86_7_n';'HD86_7_d';...
%     'TEL88_9_n';'TEL88_9_d';'HD90_1_n';'HD90_1_d';'HD92_3_n';'HD92_3_d'}; %;'TEL94_5_n';'TEL94_5_d'};

% deletion: del_amp_flag=1
% amplification: del_amp_flag=2

del_thresh = handles.com_aber_params.del_thresh; %1.6;
amp_thresh = handles.com_aber_params.amp_thresh; %2.5;
Q = handles.com_aber_params.fdr;  %0.05;
ucsc_webpage = handles.com_aber_params.ucsc_webpage;
ucsc_mounted_folder = handles.com_aber_params.ucsc_mounted_folder;

if ~exist(fullfile(workPath,'output'),'dir')
    mkdir(fullfile(workPath,'output'));
end
if isempty(ucsc_mounted_folder)
    ucsc_mounted_folder=fullfile(workPath,'output','ucsc');
    if ~exist(fullfile(workPath,'output','ucsc'),'dir')
        mkdir(fullfile(workPath,'output','ucsc'));
    end
end


params_mat={'deletion';'amplification';'del_thresh';'amp_thresh';'FDR';'ucsc_webpage';'ucsc_mounted_folder'};
params_mat(:,2)=struct2cell(handles.com_aber_params);

new_file='com_aber_params.txt';
try
    saveCellFile(params_mat,fullfile(workPath,new_file));
catch
    errordlg(['Cannot save file ' new_file ' . File may be already opened.'],'File Error','modal');
    return;
end

%if handles.com_aber_params.del==1
    del_amp_flag = 1;
    [segments_del mat_genes_del segments_samples_del sample_names_d chr_loc_and_num_samples_del]=find_common_aberr(samples, [workPath '\'], ...
        handles.ChipName, snp_ids, chr_loc_vec, chr_vec, snp_gene_symbols, snp_gene_descr, del_amp_flag, ...
        del_thresh, amp_thresh, Q, ucsc_webpage, genes_db_struct);
%end

%if handles.com_aber_params.amp==1
    del_amp_flag = 2;
    [segments_amp mat_genes_amp segments_samples_amp dum chr_loc_and_num_samples_amp]=find_common_aberr(samples, [workPath '\'], handles.ChipName, snp_ids, chr_loc_vec, ...
        chr_vec, snp_gene_symbols, snp_gene_descr, del_amp_flag, ...
        del_thresh, amp_thresh, Q, ucsc_webpage, genes_db_struct);
%end

segments=[segments_del; segments_amp];

DispStruct.SampleName='All Samples';
DispStruct.ChipName=handles.ChipName;
DispStruct.GenomeBuild=genome_assembly;
for i=1:24
    DispStruct.Chrom{i}=[];
end
chr_num=unique(segments(:,4));
for i=1:length(chr_num)
    DispStruct.Chrom{chr_num(i)}.Segments=segments(segments(:,4)==chr_num(i),1:3);
    
    %num samples del/amp
    loc_and_num_samples_del=chr_loc_and_num_samples_del(chr_loc_and_num_samples_del(:,1)==chr_num(i),2:end);
    loc_and_num_samples_amp=chr_loc_and_num_samples_amp(chr_loc_and_num_samples_amp(:,1)==chr_num(i),2:end);
    union_Locs=union(loc_and_num_samples_del(:,1),loc_and_num_samples_amp(:,1));
    DispStruct.Chrom{chr_num(i)}.Locs=union_Locs;
    DispStruct.Chrom{chr_num(i)}.Data=zeros(length(union_Locs),2);
    
    [dum idx]=ismember(loc_and_num_samples_del(:,1),union_Locs);
    DispStruct.Chrom{chr_num(i)}.Data(idx,1)=loc_and_num_samples_del(:,2);
    
    [dum idx]=ismember(loc_and_num_samples_amp(:,1),union_Locs);
    DispStruct.Chrom{chr_num(i)}.Data(idx,2)=loc_and_num_samples_amp(:,2);
    
end

save(fullfile(workPath,'output','stretches_disp.mat'),'DispStruct','params_mat');
save(fullfile(workPath,'output','genes_del.mat'),'mat_genes_del');
save(fullfile(workPath,'output','genes_amp.mat'),'mat_genes_amp');

%save gene files
table_title={'Gene' 'Chr' 'Band' 'Location' 'P-value' 'Normal Aberr.' 'UCSC' 'Num samples involved' 'Stretches Length'};
new_file=['DEL_genes_' handles.ChipName '.txt'];
try
    saveCellFile([table_title; mat_genes_del],fullfile(workPath,'output',new_file));
catch
    errordlg(['Cannot save file ' new_file ' . File may be already opened.'],'File Error','modal');
    return;
end
new_file=['AMP_genes_' handles.ChipName '.txt'];
try
    saveCellFile([table_title; mat_genes_amp],fullfile(workPath,'output',new_file));
catch
    errordlg(['Cannot save file ' new_file ' . File may be already opened.'],'File Error','modal');
    return;
end



err=make_ucsc_files(ucsc_mounted_folder, segments_del, mat_genes_del, segments_samples_del, ...
    segments_amp, mat_genes_amp, segments_samples_amp, genes_db_struct, handles.chip_annot, sample_names_d);

if ~isempty(err)
    errordlg(err,'Error in make_ucsc_files','modal');
    return;
end

msgbox('Find Common Aberrations has finished.','Find Common Aberrations Finished','modal');



% --- Executes on button press in checkbox_all_samples.
function checkbox_all_samples_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_all_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all_samples




% --------------------------------------------------------------------
function tools_item_del_genes_Callback(hObject, eventdata, handles)
% hObject    handle to tools_item_del_genes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

file_name=fullfile(workPath,'output','genes_del.mat');
if ~exist(file_name,'file')
    errordlg(['File: ' file_name ' not found'],'File Error','modal');
    return;
end

if ~isfield(handles,'chip_annot')
    errordlg('Please load data','Data not Loaded','modal');
    return;
end

load(file_name);

[selection,ok] = listdlg('ListString',cellstr([char(mat_genes_del(:,1)) char(mat_genes_del(:,3))]),'SelectionMode','single','Name', ...
    'Del. Genes', 'PromptString','Please select a gene to zoom in');


if ~ok
    return;
end

%handles = guidata(hObject);
geneSymbol = mat_genes_del{selection,1};

snpIx=strmatch(lower(geneSymbol), lower(handles.chip_annot.snps_gene_symbols), 'exact');
geneIx=strmatch(lower(geneSymbol), lower(handles.annot.symbols), 'exact');
if isempty(snpIx) || isempty(geneIx)
    errordlg(['Cannot find gene symbol: ' geneSymbol],'Gene not Found','modal');
    return;
end
xmin=min(handles.chip_annot.chr_loc_vec(min(snpIx))-100000, handles.annot.start(geneIx));
xmax=max(handles.chip_annot.chr_loc_vec(max(snpIx))+100000, handles.annot.end(geneIx));
chrom=handles.chip_annot.chr_vec(snpIx(1));

handles.h_lines_standAlone=drawChrom(handles, chrom, 1, 0);
zoomIntoChrom(handles,xmin,xmax,chrom);
guidata(gcbo, handles);





% --------------------------------------------------------------------
function tools_item_amp_genes_Callback(hObject, eventdata, handles)
% hObject    handle to tools_item_amp_genes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

file_name=fullfile(workPath,'output','genes_amp.mat');
if ~exist(file_name,'file')
    errordlg(['File: ' file_name ' not found'],'File Error','modal');
    return;
end

if ~isfield(handles,'chip_annot')
    errordlg('Please load data','Data not Loaded','modal');
    return;
end

load(file_name);

[selection,ok] = listdlg('ListString',cellstr([char(mat_genes_amp(:,1)) char(mat_genes_amp(:,3))]),'SelectionMode','single','Name', ...
    'Amp. Genes', 'PromptString','Please select a gene to zoom in');

if ~ok
    return;
end

%handles = guidata(hObject);
geneSymbol = mat_genes_amp{selection,1};

snpIx=strmatch(lower(geneSymbol), lower(handles.chip_annot.snps_gene_symbols), 'exact');
geneIx=strmatch(lower(geneSymbol), lower(handles.annot.symbols), 'exact');
if isempty(snpIx) || isempty(geneIx)
    errordlg(['Cannot find gene symbol: ' geneSymbol],'Gene not Found','modal');
    return;
end
xmin=min(handles.chip_annot.chr_loc_vec(min(snpIx))-100000, handles.annot.start(geneIx));
xmax=max(handles.chip_annot.chr_loc_vec(max(snpIx))+100000, handles.annot.end(geneIx));
chrom=handles.chip_annot.chr_vec(snpIx(1));

handles.h_lines_standAlone=drawChrom(handles, chrom, 1, 0);
zoomIntoChrom(handles,xmin,xmax,chrom);
guidata(gcbo, handles);




% --- Executes on button press in pushbutton_pca.
function pushbutton_pca_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end


chip_list=get(handles.popupmenuChipName,'string');
chip_name=chip_list{get(handles.popupmenuChipName,'value')};
if strcmp(chip_name,'Hind50k ( from Affy 100K)')
    chip_name='Hind';
elseif strcmp(chip_name,'Xba50k ( from Affy 100K)')
    chip_name='Xba';
else
    errordlg('Invalid chip','Chip Error','modal');
    return;
end

file_name=fullfile(workPath,['after_run_hmm_' chip_name '.mat']);
if ~exist(file_name,'file')
    errordlg(['Run HMM has not been performed. ' file_name ' not found.'],'Results Error','modal');
    return;
end
load(file_name);

chrom=get(handles.popupmenu_chr,'value');
if chrom==1
    chrom=-1;
else
    chrom=chrom-1;
end

DisplayAllSamplesPCA(workPath, get(handles.checkbox_genotypes,'value'), samples, chrom, 1000);
 

% --- Executes on button press in checkbox_genotypes.
function checkbox_genotypes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_genotypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_genotypes




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


% --------------------------------------------------------------------
function tools_item_display_raw_data_Callback(hObject, eventdata, handles)
% hObject    handle to tools_item_display_raw_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end
samples_from_list=get(handles.popupmenu_samples,'string');
if ischar(samples_from_list)
    errordlg('Samples should be loaded first','Samples Error','modal');
    return;
end
samples_from_list(1)=[];

chip_list=get(handles.popupmenuChipName,'string');
chip_name=get_chip_name(chip_list{get(handles.popupmenuChipName,'value')});

options=display_raw_data_GUI;
delete(options.figure1);

if options.ok==1
    if strcmp(options.params.region,'current')
        xlim1=get(handles.chromAxes,'xlim');
        region_start=round(xlim1(1));
        region_end=round(xlim1(2));
        chr_num=get(handles.chromAxes, 'UserData');
        if isempty(chr_num)
            errordlg('A chromosome should be chosen','Error','modal');
            return;
        end
            
    else
        region_start=options.params.start;
        region_end=options.params.end;
        chr_num=options.params.chrom;
    end
    
    [err, raw_copy_num_mat, chr_num_snps_vec, genotype_mat, data_snp_ids] = load_hmm_out_copy_num_mat(samples_from_list, workPath, chip_name);
    
    if ~isempty(err)
        errordlg(err,'Error','modal');
        return;
    end
    
    handles=get_chip_annot(handles,chip_name);
    
    idx=find(handles.chip_annot.chr_vec==chr_num & handles.chip_annot.chr_loc_vec>=region_start & handles.chip_annot.chr_loc_vec<=region_end);
    loc_to_plot=handles.chip_annot.chr_loc_vec(idx);
    SNPs_to_plot=handles.chip_annot.snp_ids(idx);
    [dum idx]=ismember(SNPs_to_plot, data_snp_ids);
    idx=idx(idx>0);
    raw_mat_to_plot = raw_copy_num_mat(idx,:)';
    genotype_mat_to_plot = genotype_mat(idx,:);
    [sample_pairs_cell, pairs_samples_ind] = samples_into_pairs(samples_from_list);
    genotype_mat_to_plot = geno_hmm_into_affy(genotype_mat_to_plot, SNPs_to_plot, chip_name);
    call_type_vec=zeros(size(sample_pairs_cell, 1),size(genotype_mat_to_plot,1));
    for i = 1:size(sample_pairs_cell, 1)
        call_type_vec(i,:) = call_vecs_into_call_types(genotype_mat_to_plot(:, pairs_samples_ind(i,1)), ...
            genotype_mat_to_plot(:, pairs_samples_ind(i,2)));

    end
    [loh_i loh_j]=find(call_type_vec==6);
    
    if options.params.disease_only
        raw_mat_to_plot=raw_mat_to_plot(pairs_samples_ind(:,2),:);
        samples_from_list=samples_from_list(pairs_samples_ind(:,2));
        LOH_y_to_plot=loh_i;
    else
        LOH_y_to_plot=pairs_samples_ind(loh_i,2);
    end

    
    figure, imagesc(raw_mat_to_plot);
    colorbar;
    set(gca,'xtick',[]); %1:length(SNPs_to_plot),'XTickLabel','');
    ylim1=get(gca,'ylim');
    text(1:length(SNPs_to_plot),repmat(ylim1(2)+0.01*(ylim1(2)-ylim1(1)), length(SNPs_to_plot),1),num2str(loc_to_plot),'rotation',270);
    
    set(gca,'ytick',1:length(samples_from_list));
    set(gca,'yticklabel',samples_from_list);
    
    if chr_num==23
        chr_str='X';
    elseif chr_num==24
        chr_str='Y';
    else
        chr_str=num2str(chr_num)
    end
    title(['Chr' chr_str ':' num2str(region_start) '-' num2str(region_end)]);
    
    
    text(loh_j,LOH_y_to_plot,'L');

    
    guidata(hObject, handles);
end




