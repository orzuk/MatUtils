% Prepare the .txt file for all the samples for the Hapmap SNP data.
% Here we do the Hind and Xba chips

mat=cell(271,3);

f = dir('\\eytan14-pc\Or\hapmap\data\*.CEL*');

for i=1:length(f)
    file_names{i} = f(i).name;
end

HIND_Inds = zeros(1,length(f));
XBA_Inds = zeros(1,length(f));
CEU_Inds = zeros(1,length(f));
JPT_CHB_Inds = zeros(1,length(f));
YRI_Inds = zeros(1,length(f));

for i=1:length(f)
    if(strfind(file_names{i}, 'HIND'))
        HIND_Inds(i) = 1;
    end
    if(strfind(file_names{i}, 'XBA'))
        XBA_Inds(i) = 1;
    end
end
HIND_Inds = find(HIND_Inds);
XBA_Inds = find(XBA_Inds);
CEU_Inds = strmatch('CEU', file_names);
JPT_CHB_Inds = union( strmatch('JPT', file_names), strmatch('CHB', file_names));
YRI_Inds = strmatch('YRI', file_names);


load files.mat;
load AllHapmapGenders.mat;

% Make the JPT always before CHB
for i=1:length(file_names)
    correct_file_names{i} = file_names{i};
end
for i=JPT_CHB_Inds'
    correct_file_names{i}(1:7) = 'JPT_CHB';
end

for i=JPT_CHB_Inds'
    i_is = i
    sys_str = ['''move ' fullfile('\\eytan14-pc\Or\hapmap\data\', file_names{i})  '  '  ...
        fullfile('\\eytan14-pc\Or\hapmap\data\', correct_file_names{i})  '''' ];
    eval(['system(' sys_str ');']);
end


for i=1:length(correct_file_names)
    idx=strfind(correct_file_names{i},'_NA');
    correct_sample_names{i} = correct_file_names{i}(1:idx+7);
    
    if(strfind(correct_sample_names{i}, 'NA18996'))
        correct_sample_names{i}(end-6:end) = 'NA19012';
    end
    
end
    
[SAMPLE_INTERSECT III JJJ] = intersect(correct_sample_names, AllHapMapSamples.Name)

% Difference:  JPT_CHB_NA18996 only in correct_sample_names (in Affy Chips)
% Difference:  JPT_CHB_NA19012 only in AllHapMapSamples.Name (in Genotypes data) probably the correct one


%hind
mat{1,1}='sample';
mat{1,2}='CEL file';
mat{1,3}='Gender';

for i=1:length(HIND_Inds)
    i
    idx=strfind(correct_file_names{HIND_Inds(i)},'_HIND');
    mat{i+1,1}=[correct_file_names{HIND_Inds(i)}(1:idx-1) '_n'];
    mat{i+1,2}=fullfile('f:\Or\hapmap\data',correct_file_names{HIND_Inds(i)});    
    GenderIdx = strmatch(correct_sample_names{HIND_Inds(i)}, AllHapMapSamples.Name);
    mat{i+1,3}=AllHapMapSamples.Gender{GenderIdx}; %  'M';
end
saveCellFile(mat,'samples_hind.txt');

%xba
for i=1:length(XBA_Inds)
    i
    idx=strfind(correct_file_names{XBA_Inds(i)},'_XBA');
    mat{i+1,1}=[correct_file_names{XBA_Inds(i)}(1:idx-1) '_n'];
    mat{i+1,2}=fullfile('f:\Or\hapmap\data',correct_file_names{XBA_Inds(i)});    
    GenderIdx = strmatch(correct_sample_names{XBA_Inds(i)}, AllHapMapSamples.Name);
    mat{i+1,3}=AllHapMapSamples.Gender{GenderIdx}; %  'M';
end
saveCellFile(mat,'samples_xba.txt');
