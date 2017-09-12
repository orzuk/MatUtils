function handles=get_chip_annot(handles,chip_name)

if strcmpi(chip_name,'Hind')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Hind')
        handles.chip_annot=load(fullfile('..','database',['Hind_annot_data_' get_genome_assembly() '.mat']));
        handles.ChipName='Hind';
    end
elseif strcmpi(chip_name,'Xba')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Xba')
        handles.chip_annot=load(fullfile('..','database',['Xba_annot_data_' get_genome_assembly() '.mat']));
        handles.ChipName='Xba';
    end
elseif strcmp(chip_name,'Hind_and_Xba')
    if ~isfield(handles,'ChipName') || ~strcmpi(handles.ChipName,'Hind_and_Xba')
        handles.chip_annot=load_Hind_and_Xba;
        handles.ChipName='Hind_and_Xba';
    end
    % else
    %     errordlg(['Invalid chip ' chip_name],'Chip Error','modal');
    %     return;
end