% Prepare the .txt file for all the samples for the Colon SNP data.
% Here everything is Xba.

f = dir('E:\Or\colon\*.CEL');

mat=cell(length(f)+1,3);

for i=1:length(f)
    file_names{i} = f(i).name;
end

% From here we'll get the geneders
T = loadCellFile('E:\Or\colon\SNP samples details_061906_PLUS_CEL.txt');
T_labels = T(:,1); T_file_names = T(:,5); T_gender = T(:,8);

%xba
mat{1,1}='sample';
mat{1,2}='CEL file';
mat{1,3}='Gender';

EmptyIdxVec=zeros(1,length(file_names));
file_tissue_names = {};
for i=1:length(file_names)
    i
    idx=strfind(file_names{i},'.CEL');
    sample_type = file_names{i}(end-4);
    if(sample_type == 'H') % Healthy individual
        mat{i+1,1}=[file_names{i}(1:idx-1) '_n'];
    else
        mat{i+1,1}=[file_names{i}(1:idx-1) '_d'];
    end
    mat{i+1,2}=fullfile('E:\Or\colon',file_names{i});

    GenderIdx = strmatch(file_names{i},T_file_names);
    mat{i+1,3}=T_gender{GenderIdx}; 
end

saveCellFile(mat,'colon_samples_xba.txt');

