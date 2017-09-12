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

% Last Modified by GUIDE v2.5 30-Jul-2007 16:04:39

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
AssignAllGlobalConstants(); %% genome_assembly = get_genome_assembly();
load(fullfile('..','database',['chr_data_' genome_assembly '.mat']));
handles.chr_band_end = chr_band_end;
handles.chr_band_start = chr_band_start;
handles.chr_band_names = chr_band_names;
handles.chr_num = chr_num;
handles.end_p_location = end_p_location;
load(fullfile('..','database',['refgenes_' genome_assembly '.mat']));
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
handles.com_aber_params_default.del_thresh=1.85;
handles.com_aber_params_default.amp_thresh=2.5;
handles.com_aber_params_default.fdr=0.05;
handles.com_aber_params_default.min_stretch=3;
handles.com_aber_params_default.max_stretch=10;
handles.com_aber_params_default.disease=1;
handles.com_aber_params_default.normals=0;
handles.com_aber_params_default.remove_del_arm=0;
handles.com_aber_params_default.remove_amp_arm=0;
handles.com_aber_params_default.fraction_arm=0.3;
%handles.com_aber_params_default.ucsc_webpage='http://www.weizmann.ac.il/';
handles.com_aber_params_default.ucsc_webpage='http://www.weizmann.ac.il/complex/Domany/libi/Leukemia_all_ucsc';
handles.com_aber_params_default.ucsc_mounted_folder='';

handles.com_aber_params=handles.com_aber_params_default;

set(handles.popupmenu_samples,'userdata',get(handles.popupmenu_samples,'string'));


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
chip_name=get_chip_name(chip_name);

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

mat=loadCellFile_str(samplesFile);
if size(mat,2)<3 || ~strcmpi(mat{1,1},'sample') || ~strcmpi(mat{1,2},'cel file')  || ~strcmpi(mat{1,3},'gender')
    errordlg({'First column title must be "Sample"' 'Second column title must be "CEL file"' 'Third column title must be "Gender"'},'File Error','modal');
    return;
end

if size(mat,2)==4 && ~strcmpi(mat{1,4},'diploid')
    errordlg('The optional fourth column title must be "Diploid"','File Error','modal');
    return;
end
    
samples=mat(2:end,1);
arrays=mat(2:end,2);
gender=mat(2:end,3);
gender_not_given_vec=find(strcmpi(gender,'nan'));
if length( find(strcmp(gender,'M') | strcmp(gender,'F') | strcmpi(gender,'nan')) ) ~= length(gender)
    errordlg('The 3rd column values in Samples file must be either M, F or NAN','Input Error','modal');
    return;
end

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
    answer=inputdlg({'Number of EM iterations' 'Number of EM starting points' 'EM_tolerance' 'Genotype correlations (0,1,2)'} ...
        ,'HMM Par.',1,...
        {'25' '5' '0.000001' '2'});

    if isempty(answer)
        return;
    end
    handles.HMMParamsStruct.num_EM_iters = str2num(answer{1}); %50;
    handles.HMMParamsStruct.num_EM_starting_points = str2num(answer{2}); %10;
    handles.HMMParamsStruct.EM_tolerance = str2num(answer{3}); %0.000001;
    switch (answer{4})
        case '0'
            handles.HMMParamsStruct.UseGenotypesCorrs=0;
        case '1'
            handles.HMMParamsStruct.UseGenotypesCorrs=1;
        case '2'
            handles.HMMParamsStruct.UseGenotypesCorrs=2;
        otherwise
            errordlg('"Use genotype correlations" field must be either 0,1 or 2','Parameter Error','modal');
            return;
    end
    handles.run_HMM_after_load=1;
else
    handles.run_HMM_after_load=0;
end


if size(mat,2)==4
    normal_ind = find(strcmp(mat(2:end,4),'2'));
else
    normal_ind = get_normal_copy_num_samples_ind(samples);
end
female_ind = strmatch('F', gender); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_cel_files_flag=1; prep_block_flag=1; % Can change to zero after first time ..
female_ind=RLMM_Normalize(workPath, samples, normal_ind,  female_ind, arrays, chip_name, load_cel_files_flag, ...
    prep_block_flag, gender_not_given_vec); % Do not give correct dir (it will overide files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gender(:)=cellstr('M');
gender(female_ind)=cellstr('F');

diploid_ind = get_diploid_ind(samples);
%load_and_norm_data(samples, arrays, workPath, diploid_ind, chip_name);

save(fullfile(workPath,['after_normalization_' chip_name '.mat']),'samples','gender');

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

AssignAllGlobalConstants(); %% genome_assembly = get_genome_assembly();
chip_list=get(handles.popupmenuChipName,'string');
chip_name=chip_list{get(handles.popupmenuChipName,'value')};

workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

% samplesFile=deblank(get(handles.edit_samp_file,'string'));
% if exist(samplesFile)~=2
%     errordlg('Samples file not found','File Error','modal');
%     return;
% end
% 
% mat=loadCellFile_str(samplesFile);
% if size(mat,2)<3 || ~strcmpi(mat{1,1},'sample') || ~strcmpi(mat{1,2},'cel file')  || ~strcmpi(mat{1,3},'gender')
%     errordlg({'First column title must be "Sample"' 'Second column title must be "CEL file"' 'Third column title must be "Gender"'},'File Error','modal');
%     return;
% end
% samples=mat(2:end,1);
% gender=mat(2:end,3);
% if length(find(strcmpi(gender,'M') | strcmpi(gender,'F')))~= length(gender)
%     errordlg('The 3rd column values in Samples file must be either M or F','Input Error','modal');
%     return;
% end

if isfield(handles,'run_HMM_after_load') && handles.run_HMM_after_load
    HMMParamsStruct=handles.HMMParamsStruct;
else
    answer=inputdlg({'Number of EM iterations' 'Number of EM starting points' 'EM_tolerance' 'Genotype correlations (0,1,2)'} ...
        ,'HMM Par.',1,...
        {'25' '5' '0.000001' '2'});

    if isempty(answer)
        return;
    end
    HMMParamsStruct.num_EM_iters = str2num(answer{1}); %50;
    HMMParamsStruct.num_EM_starting_points = str2num(answer{2}); %10;
    HMMParamsStruct.EM_tolerance = str2num(answer{3}); %0.000001;

    switch (answer{4})
        case '0'
            HMMParamsStruct.UseGenotypesCorrs=0;
        case '1'
            HMMParamsStruct.UseGenotypesCorrs=1;
        case '2'
            HMMParamsStruct.UseGenotypesCorrs=2;
        otherwise
            errordlg('"Use genotype correlations" field must be either 0,1 or 2','Parameter Error','modal');
            return;
    end    
end
% what is hapmap_population ?
if strcmp(chip_name,'Hind50k ( from Affy 100K)')
    LDStruct=load(fullfile('..','database',['LD_' pop_str_vec{default_hapmap_population} '_hind_' genome_assembly '.mat']));
    handles.chip_annot=load(fullfile('..','database',['Hind_annot_data_' genome_assembly '.mat']));
    handles.ChipName='Hind';
elseif strcmp(chip_name,'Xba50k ( from Affy 100K)')
    LDStruct=load(fullfile('..','database',['LD_' pop_str_vec{default_hapmap_population} '_xba_' genome_assembly '.mat']));
    handles.chip_annot=load(fullfile('..','database',['Xba_annot_data_' genome_assembly '.mat']));
    handles.ChipName='Xba';
else
    errordlg('Invalid chip','Chip Error','modal');
    return;
end

if(HMMParamsStruct.num_EM_iters == -1) % special case: RLMM learning
    UpdateDatabaseFlag=1;
    workPath2 = workPath; % fullfile(workPath, 'CEU_test');
else
    UpdateDatabaseFlag=0;
    workPath2 = workPath;
end

file_name=fullfile(workPath2,['after_normalization_' handles.ChipName '.mat']);
if ~exist(file_name,'file')
    errordlg(['Normalization has not been performed. ' file_name ' not found.'],'Results Error','modal');
    return;
end
load(file_name);
    
HMMParamsStruct.x_dim = 3; HMMParamsStruct.y_dim = 1;


HMMParamsStruct.learn_model_EM = 1;
HMMParamsStruct.ChromosomesToRun = 1:num_chr_to_run;
HMMParamsStruct.ChromosomesToPlot = 1:num_chr_to_run; % Plot

HMMParamsStruct.derich = 1.0/120.0;

    
%% OutputFilesNames = FakeAfterRunHMMFile(workPath, samples, LDStruct, handles.chip_annot, HMMParamsStruct);

HMMParamsStruct.hapmap_dir = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data';
% for pop=1:3
%     if(~isempty(findstr(pop_str_vec{pop}, samplesFile)))
%         HMMParamsStruct.hapmap_population = pop;
%     end
% end
HMMParamsStruct.hapmap_version = genome_assembly_to_hapmap_version(); % 'r22'
HMMParamsStruct.hapmap_population = CEU; 

if(UpdateDatabaseFlag)
    MinSparse=5; PsuedoCount=40; 
    HMMParamsStruct.hapmap_population = CEU; % ALL_POPS; % Force to take all populations

    if(HMMParamsStruct.hapmap_population == ALL_POPS)
        for hapmap_pop2 = 1:3   
            workPath2 = fullfile(workPath, [pop_str_vec{hapmap_pop2} '_test']);
            [OutputFilesNames NormalMeanIntensities r_mat] = ...
                PreLearnMoGChromFromData(workPath2, lower(handles.chip_annot.chip));
        end
    end

    OutputFilesName = UpdateDatabaseRLMMParams(workPath, samples, LDStruct, handles.chip_annot, ...
        HMMParamsStruct, HMMParamsStruct.hapmap_dir, HMMParamsStruct.hapmap_population, HMMParamsStruct.hapmap_version, ...
        MinSparse, PsuedoCount);
else
    [OutputFilesNames NormalMeanIntensities r_mat] = ...
        PreLearnMoGChromFromData(workPath, lower(handles.chip_annot.chip));

    HMMParamsStruct.r_mat = r_mat; % copy the r relative parameters
    OutputFilesNames = LearnHMMChromFromData(workPath, samples, LDStruct, handles.chip_annot, HMMParamsStruct);
    % handles.OutputFilesNames=OutputFilesNames;
    guidata(hObject, handles);

    save(fullfile(workPath,['after_run_hmm_' handles.ChipName '.mat']),'OutputFilesNames','samples','gender');
end

set(handles.popupmenu_samples,'string',['Average'; get_samples_with_gender(samples, gender)]);
set(handles.popupmenu_samples,'userdata',['Average'; samples]);

msgbox('Run HMM has finished.','Run HMM Finished','modal');


% --- Executes on button press in pushbutton_view_res.
function pushbutton_view_res_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_view_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%genome_assembly = get_genome_assembly();
workPath = deblank(get(handles.editWorkPath,'string'));
if ~isdir(workPath)
    errordlg('Working path not found','Path Error','modal');
    return;
end

chip_list=get(handles.popupmenuChipName,'string');
chip_name=chip_list{get(handles.popupmenuChipName,'value')};
chip_name=get_chip_name(chip_name);
PairedSample=[];

if get(handles.checkbox_all_samples','value')
    if handles.com_aber_params.run_chips_together
        chip_name='Hind_and_Xba';
    end
    if handles.com_aber_params.disease
        normal_disease_suffix='_d';
    else
        normal_disease_suffix='_n';
    end
    file_name=['stretches_disp_' chip_name normal_disease_suffix '.mat'];
    if ~exist(fullfile(workPath,'output',file_name),'file')
        errordlg(['Find Common Aberrations has not been performed. ' file_name ' not found.'],'Results Error','modal');
        return;
    end
    load(fullfile(workPath,'output',file_name));
else
    
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

    samples_from_list=get(handles.popupmenu_samples,'userdata');
    samples_from_list(1)=[];
    if ~isequal(samples_from_list,samples)
        set(handles.popupmenu_samples,'string',['Average'; get_samples_with_gender(samples, gender)]);
        set(handles.popupmenu_samples,'userdata',['Average'; samples]);
        msgbox('Please choose a sample and click "View Results" again.','View Results','modal');
        return;
    end
    
    %load(OutputFilesNames.disp_files{get(handles.popupmenu_samples,'value')});
    load(fullfile(workPath,OutputFilesNames.disp_files{get(handles.popupmenu_samples,'value')}));
    
    %if there is a paired sample
    if strcmpi(DispStruct.SampleName(end-1:end),'_d')
        idx=find(strcmpi([DispStruct.SampleName(1:end-2) '_n'],samples));
        if ~isempty(idx)
            PairedSample=load(fullfile(workPath,OutputFilesNames.disp_files{idx+1}));
        end
    end
end

handles=get_chip_annot(handles,DispStruct.ChipName);

%load data\SampleChromsDataStructForDisplay.mat;
handles.GenomeBuild=DispStruct.GenomeBuild;
handles.SampleName=DispStruct.SampleName;
handles.ChromsData=DispStruct.Chrom;

handles.PairedSample=PairedSample;
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
if ~isempty(find(~cellfun('isempty',handles.ChromsData)))
    [handles.subplots handles.h_lines handles.min_intensity_seg_all handles.max_intensity_seg_all]=drawChromMap(handles, 0);
else
    msgbox('Data are empty.','View Results','modal');
    return;
end
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
xmin1=min(handles.chip_annot.chr_loc_vec(min(snpIx))-100000, handles.annot.start(geneIx));
xmax1=max(handles.chip_annot.chr_loc_vec(max(snpIx))+100000, handles.annot.end(geneIx));

xmin=xmin1-(xmax1-xmin1);
xmax=xmax1+(xmax1-xmin1);

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

mat=loadCellFile_str(samplesFile);
if size(mat,2)<3 || ~strcmpi(mat{1,1},'sample') || ~strcmpi(mat{1,2},'cel file')  || ~strcmpi(mat{1,3},'gender')
    errordlg({'First column title must be "Sample"' 'Second column title must be "CEL file"' 'Third column title must be "Gender"'},'File Error','modal');
    return;
end

samples=mat(2:end,1);
gender=mat(2:end,3);
gender_not_given_vec=find(strcmpi(gender,'nan'));
if length( find(strcmp(gender,'M') | strcmp(gender,'F') | strcmpi(gender,'nan')) ) ~= length(gender)
    errordlg('The 3rd column values in Samples file must be either M, F or NAN','Input Error','modal');
    return;
end

if handles.com_aber_params.run_chips_together
    chip_name='Hind_and_Xba';
end

chip_name=get_chip_name(chip_name);
handles=get_chip_annot(handles,chip_name);

guidata(hObject, handles);

if ~isempty(gender_not_given_vec)
    after_run_hmm=load(fullfile(workPath,['after_run_hmm_' chip_name '.mat']));
    [dum idx]=ismember(samples(gender_not_given_vec),after_run_hmm.samples)
    gender(gender_not_given_vec)=after_run_hmm.gender(idx);
end
    



AssignAllGlobalConstants();  

genes_db_f = ['refgenes_' genome_assembly '.mat'];

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
min_stretch = handles.com_aber_params.min_stretch;
max_stretch = handles.com_aber_params.max_stretch;
disease = handles.com_aber_params.disease;
normals = handles.com_aber_params.normals;
remove_del_arm = handles.com_aber_params.remove_del_arm;
remove_amp_arm = handles.com_aber_params.remove_amp_arm;
fraction_arm = handles.com_aber_params.fraction_arm;
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


params_mat={'deletion';'amplification';'del_thresh';'amp_thresh';'FDR';'min_stretch';'max_stretch';'disease';'normals';...
    'remove_del_arm';'remove_amp_arm';'fraction_arm';'ucsc_webpage';'ucsc_mounted_folder'};
params_mat(:,2)=struct2cell(handles.com_aber_params);

new_file='com_aber_params.txt';
try
    saveCellFile(params_mat,fullfile(workPath,new_file));
catch
    errordlg(['Cannot save file ' new_file ' . File may be already opened.'],'File Error','modal');
    return;
end


hmm_out=load_hmm_out_prepare_zero_one_mat_of_chips(samples, workPath, handles.ChipName, del_thresh, amp_thresh, ...
    snp_ids, handles.chip_annot.chr_vec, gender);

normal_disease_flag_vec=[];
normal_disease_suffix_vec={};
if disease
    normal_disease_flag_vec(end+1)=2;
    normal_disease_suffix_vec{end+1}='_d';
end
if normals
    normal_disease_flag_vec(end+1)=1;
    normal_disease_suffix_vec{end+1}='_n';
end
%remove_outlier_flag = 1;
remove_outlier_flag = 0;
outlier_thresh = 0.1;
if(remove_outlier_flag)
    [hmm_out, good_ind] = remove_outliers_hmm_out(hmm_out, outlier_thresh);
    samples = samples(good_ind);
    gender = gender(good_ind);
end
    
for j=1:length(normal_disease_flag_vec)

    del_amp_flag = 1;
    [segments_del mat_genes_del segments_samples_del sample_names_d chr_loc_and_num_samples_del segments_samples_arm_del]=...
        find_common_aberr(samples, gender, workPath, handles, snp_ids, chr_loc_vec, chr_vec, ...
        snp_gene_symbols, snp_gene_descr, del_amp_flag, Q, min_stretch, max_stretch, normal_disease_flag_vec(j), remove_del_arm, ...
        fraction_arm, ucsc_webpage, genes_db_struct, gene_symbols, gene_descr, hmm_out);

    del_amp_flag = 2;
    [segments_amp mat_genes_amp segments_samples_amp dum chr_loc_and_num_samples_amp segments_samples_arm_amp]=...
        find_common_aberr(samples, gender, workPath, handles, snp_ids, chr_loc_vec, chr_vec, ...
        snp_gene_symbols, snp_gene_descr, del_amp_flag, Q, min_stretch, max_stretch, normal_disease_flag_vec(j), remove_amp_arm, ...
        fraction_arm, ucsc_webpage, genes_db_struct, gene_symbols, gene_descr, hmm_out);


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
        
        %add segments edges with values zero
        Locs_to_add=DispStruct.Chrom{chr_num(i)}.Segments(:,1)-1;
        Locs_to_add=[Locs_to_add; DispStruct.Chrom{chr_num(i)}.Segments(:,2)+1];
        
        DispStruct.Chrom{chr_num(i)}.Locs = [DispStruct.Chrom{chr_num(i)}.Locs; Locs_to_add];
        DispStruct.Chrom{chr_num(i)}.Data = [DispStruct.Chrom{chr_num(i)}.Data; zeros(length(Locs_to_add),2)];
        
        [DispStruct.Chrom{chr_num(i)}.Locs idx] = sort(DispStruct.Chrom{chr_num(i)}.Locs);
        DispStruct.Chrom{chr_num(i)}.Data = DispStruct.Chrom{chr_num(i)}.Data(idx,:);
        
    end

    save(fullfile(workPath,'output',['stretches_disp_' handles.ChipName normal_disease_suffix_vec{j} '.mat']),'DispStruct','params_mat');
    save(fullfile(workPath,'output',['genes_del_' handles.ChipName normal_disease_suffix_vec{j} '.mat']),'mat_genes_del');
    save(fullfile(workPath,'output',['genes_amp_' handles.ChipName normal_disease_suffix_vec{j} '.mat']),'mat_genes_amp');

    %save gene files
    table_title={'Gene' 'Gene Descr.' 'Cancer involvement' 'Chr' 'Band' 'Location' 'P-value' 'Normal Aberr.' 'Genecards' 'UCSC' 'Stretches' 'Num samples involved' ...
        'Stretches Length'};
    new_file=fullfile(workPath,'output',['DEL_genes_' handles.ChipName normal_disease_suffix_vec{j} '.txt']);
    if ~isempty(mat_genes_del)
        [dum idx]=sort(cell2mat(mat_genes_del(:,7))); %sort by p-value 17/07/2007, 07/08/2007
        mat_genes_del=mat_genes_del(idx,:);
        try
            saveCellFile([table_title; mat_genes_del],new_file);
        catch
            errordlg(['Cannot save file ' new_file ' . File may be already opened.'],'File Error','modal');
            return;
        end
    elseif exist(new_file,'file')
        delete(new_file);
        if exist(new_file,'file')
            errordlg(['Error deleting ' new_file],'Error','modal');
            return;
        end
    end
    new_file=fullfile(workPath,'output',['AMP_genes_' handles.ChipName normal_disease_suffix_vec{j} '.txt']);
    if ~isempty(mat_genes_amp)
        [dum idx]=sort(cell2mat(mat_genes_amp(:,7))); %sort by p-value 17/07/2007, 07/08/2007
        mat_genes_amp=mat_genes_amp(idx,:);
        try
            saveCellFile([table_title; mat_genes_amp], new_file);
        catch
            errordlg(['Cannot save file ' new_file ' . File may be already opened.'],'File Error','modal');
            return;
        end
    elseif exist(new_file,'file')
        delete(new_file);
        if exist(new_file,'file')
            errordlg(['Error deleting ' new_file],'Error','modal');
            return;
        end 
    end

    err=make_ucsc_files(ucsc_mounted_folder, segments_del, mat_genes_del, segments_samples_del, segments_samples_arm_del, ...
        segments_amp, mat_genes_amp, segments_samples_amp, segments_samples_arm_amp, ...
        genes_db_struct, handles.chip_annot, sample_names_d, normal_disease_suffix_vec{j});

    if ~isempty(err)
        errordlg(err,'Error in make_ucsc_files','modal');
        return;
    end

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

if handles.com_aber_params.disease
    normal_disease_suffix='_d';
else
    normal_disease_suffix='_n';
end

file_name=fullfile(workPath,'output',['genes_del_' handles.ChipName normal_disease_suffix '.mat']); 
if ~exist(file_name,'file')
    errordlg(['File: ' file_name ' not found'],'File Error','modal');
    return;
end

if ~isfield(handles,'chip_annot')
    errordlg('Please load data','Data not Loaded','modal');
    return;
end

load(file_name);

[selection,ok] = listdlg('ListString',cellstr([char(mat_genes_del(:,1)) repmat(' ',size(mat_genes_del,1),10) char(mat_genes_del(:,5)) ...
    repmat(' ',size(mat_genes_del,1),10) num2str(cell2mat(mat_genes_del(:,7)))]),...
    'SelectionMode','single','Name', 'Del. Genes', 'PromptString','Please select a gene to zoom in','ListSize',[240 450]);


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

if handles.com_aber_params.disease
    normal_disease_suffix='_d';
else
    normal_disease_suffix='_n';
end

file_name=fullfile(workPath,'output',['genes_amp_' handles.ChipName normal_disease_suffix '.mat']);
if ~exist(file_name,'file')
    errordlg(['File: ' file_name ' not found'],'File Error','modal');
    return;
end

if ~isfield(handles,'chip_annot')
    errordlg('Please load data','Data not Loaded','modal');
    return;
end

load(file_name);

[selection,ok] = listdlg('ListString',cellstr([char(mat_genes_amp(:,1)) repmat(' ',size(mat_genes_amp,1),10) char(mat_genes_amp(:,5))...
    repmat(' ',size(mat_genes_amp,1),10) num2str(cell2mat(mat_genes_amp(:,7)))]),...
    'SelectionMode','single','Name', 'Amp. Genes', 'PromptString','Please select a gene to zoom in', 'ListSize',[240 450]);

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
chip_name=get_chip_name(chip_name);

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

DisplayAllSamplesPCA(workPath, chip_name, get(handles.checkbox_genotypes,'value'), samples, chrom, 1000);
 

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
samples_from_list=get(handles.popupmenu_samples,'userdata');
if ischar(samples_from_list)
    errordlg('Samples should be loaded first','Samples Error','modal');
    return;
end
samples_from_list(1)=[];

samplesFile=deblank(get(handles.edit_samp_file,'string'));
if ~isempty(samplesFile)
    if ~exist(samplesFile,'file')
        errordlg('Samples file not found','File Error','modal');
        return;
    end
    mat=loadCellFile_str(samplesFile);
    if size(mat,2)<1 || ~strcmpi(mat{1,1},'sample')
        errordlg('First column title must be "Sample"','File Error','modal');
        return;
    end
    samples=mat(2:end,1);
    idx=ismember(samples_from_list,samples);
    if ~any(idx)
        errordlg('Input samples do not exist','Samples Error','modal');
        return;
    end
    samples_from_list=samples_from_list(idx);
end

if handles.com_aber_params.run_chips_together
    chip_name='Hind_and_Xba';
else
    chip_list=get(handles.popupmenuChipName,'string');
    chip_name=get_chip_name(chip_list{get(handles.popupmenuChipName,'value')});
end

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
    
    normal_disease_suffix_vec={};
    i_vec=[];
    if  options.params.normal
        normal_disease_suffix_vec{end+1}='_n';
        i_vec(end+1)=1;
    end
    if options.params.disease
        normal_disease_suffix_vec{end+1}='_d';
        i_vec(end+1)=2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot stretches
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    flag_loaded=0;
    for i=1:length(normal_disease_suffix_vec)

        del_stretch_f = fullfile(workPath, 'output', ['DEL_stretches_' chip_name normal_disease_suffix_vec{i} '.xls']);
        amp_stretch_f = fullfile(workPath, 'output', ['AMP_stretches_' chip_name normal_disease_suffix_vec{i} '.xls']);
        del_stretch_cell = [];
        if(exist(del_stretch_f,'file'))
            del_stretch_cell = loadCellFile_str(del_stretch_f);
        end
        amp_stretch_cell = [];
        if(exist(amp_stretch_f,'file'))
            amp_stretch_cell= loadCellFile_str(amp_stretch_f);
        end
        if (length(amp_stretch_cell) || length(del_stretch_cell))
            if ~flag_loaded
                genome_assembly = get_genome_assembly();
                load(fullfile('..','database',['refgenes_' genome_assembly '.mat']), ...
                    'exon_start_cell', 'exon_end_cell');
                flag_loaded=1;
            end
            %libi: maybe save 'exon_start_cell', 'exon_end_cell' in handles so that
            %it won't be loaded each time separately

            plot_common_aberr(del_stretch_cell, amp_stretch_cell, chr_num, region_start, region_end, ...
                handles.annot.chr, handles.annot.start, handles.annot.symbols, exon_start_cell, ...
                exon_end_cell, handles.chip_annot.chr_vec, handles.chip_annot.chr_loc_vec);

        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END Plot stretches
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles=get_chip_annot(handles,chip_name);
    
    if ~options.params.hmm_data % don't load HMM raw copy number - it is not used.
        [err, genotype_mat, data_snp_ids, chr_num_snps_vec] = ...
            load_hmm_out_copy_num_mat_of_chips(samples_from_list, workPath, chip_name, handles.chip_annot.snp_ids);
    else
        [err, genotype_mat, data_snp_ids, chr_num_snps_vec, raw_copy_num_mat] = ...
            load_hmm_out_copy_num_mat_of_chips(samples_from_list, workPath, chip_name, handles.chip_annot.snp_ids);
    end
    
    if ~isempty(err)
        errordlg(err,'Error','modal');
        return;
    end
    
    if ~options.params.hmm_data
        [err, raw_copy_num_mat2, data_snp_ids2] = ...
            load_normalized_copy_num_mat_of_chips(samples_from_list, workPath, chip_name, handles.chip_annot.snp_ids);
        if ~isempty(err)
            errordlg(err,'Error','modal');
            return;
        end
        [dum idx]=ismember(data_snp_ids, data_snp_ids2);
        raw_copy_num_mat=raw_copy_num_mat2(idx,:);
        clear raw_copy_num_mat2 data_snp_ids2
    end
    
    idx=find(handles.chip_annot.chr_vec==chr_num & handles.chip_annot.chr_loc_vec>=region_start & handles.chip_annot.chr_loc_vec<=region_end);
    loc_to_plot=handles.chip_annot.chr_loc_vec(idx);
    SNPs_to_plot=handles.chip_annot.snp_ids(idx);
    [dum idx]=ismember(SNPs_to_plot, data_snp_ids);
    idx=idx(idx>0);
    loc_to_plot(dum==0) = [];
    SNPs_to_plot(dum==0) = [];
    raw_mat_to_plot = raw_copy_num_mat(idx,:)';
    clear raw_copy_num_mat;
    genotype_mat_to_plot = genotype_mat(idx,:);
    clear genotype_mat;
    [sample_pairs_cell, pairs_samples_ind] = samples_into_pairs(samples_from_list);
    genotype_mat_to_plot = geno_hmm_into_affy(genotype_mat_to_plot, SNPs_to_plot, chip_name, handles.chip_annot.snp_ids);
    call_type_vec=zeros(size(sample_pairs_cell, 1),size(genotype_mat_to_plot,1));
    for i = 1:size(sample_pairs_cell, 1)
        call_type_vec(i,:) = call_vecs_into_call_types(genotype_mat_to_plot(:, pairs_samples_ind(i,1)), ...
            genotype_mat_to_plot(:, pairs_samples_ind(i,2)));

    end
    [loh_i loh_j]=find(call_type_vec==6);
    
    
%     if options.params.disease && options.params.normal
%         LOH_y_to_plot=pairs_samples_ind(loh_i,2); %plot together
%     end
        
    LOH_y_to_plot=loh_i;
    for i=i_vec

        raw_mat_to_plot2=raw_mat_to_plot(pairs_samples_ind(:,i),:);
        samples_from_list2=samples_from_list(pairs_samples_ind(:,i));

        figure, imagesc(raw_mat_to_plot2,[0 5]);
        h=colorbar;
        %ylim1=get(h,'ylim');
        %set(h,'ytick',ylim1(1)+(ylim1(2)-ylim1(1))*(0:5)/5); %[ylim1(1) mean(ylim1) ylim1(2)]);
        %set(h,'YTickLabel',{'0' '1' '2' '3' '4' '5'});
        set(gca,'xtick',[]); %1:length(SNPs_to_plot),'XTickLabel','');
        ylim1=get(gca,'ylim');
        text(1:length(SNPs_to_plot),repmat(ylim1(2)+0.01*(ylim1(2)-ylim1(1)), length(SNPs_to_plot),1),num2str(loc_to_plot),'rotation',270,'fontsize',7);

        set(gca,'ytick',1:length(samples_from_list2));
        set(gca,'yticklabel',samples_from_list2,'fontsize',7);

        if chr_num==23
            chr_str='X';
        elseif chr_num==24
            chr_str='Y';
        else
            chr_str=num2str(chr_num);
        end
        title(['Chr' chr_str ':' num2str(region_start) '-' num2str(region_end)]);


        text(loh_j,LOH_y_to_plot,'L');
    end
    
    guidata(hObject, handles);
end


% --- Executes on selection change in popupmenu_zoom_out.
function popupmenu_zoom_out_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_zoom_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_zoom_out contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_zoom_out


chrom=get(handles.chromAxes,'UserData');
if isempty(chrom)
    errordlg('A chromosome segment should be chosen.','Chromosome Error','modal');
    return;
end   
xlim1=get(handles.chromAxes,'xlim');
list1=get(handles.popupmenu_zoom_out,'string');
zoom_out_fold=str2num(list1{get(handles.popupmenu_zoom_out,'value')});
increase_ratio=(zoom_out_fold-1)/2;

xmin=xlim1(1)-diff(xlim1)*increase_ratio;
xmax=xlim1(2)+diff(xlim1)*increase_ratio;

%handles.h_lines_standAlone=drawChrom(handles, chrom, 1, 0);
zoomIntoChrom(handles,xmin,xmax,chrom);
%guidata(gcbo, handles);
    


% --- Executes during object creation, after setting all properties.
function popupmenu_zoom_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_zoom_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function new_genome_build_item_Callback(hObject, eventdata, handles)
% hObject    handle to new_genome_build_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp(1);
add_genome_build_GUI;


% --------------------------------------------------------------------
function new_chip_item_Callback(hObject, eventdata, handles)
% hObject    handle to new_chip_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function update_annotations_item_Callback(hObject, eventdata, handles)
% hObject    handle to update_annotations_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function database_sub_menu_Callback(hObject, eventdata, handles)
% hObject    handle to database_sub_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
