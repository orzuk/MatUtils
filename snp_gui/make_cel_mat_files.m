function make_CEL_mat_files(CEL_folder)

files=dir(fullfile(CEL_folder, '*.CEL'));

for i=1:length(files)
    CELStruct=affyread(fullfile(CEL_folder,files(i).name));
    CELStruct.Probes=single(CELStruct.Probes);
    save(fullfile(CEL_folder,[files(i).name(1:end-4) '.mat']),'CELStruct');
    clear CELStruct
end