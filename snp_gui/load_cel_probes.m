% Read a Cel files and output it's probs name, snp_ids 
% and the data (intensities) itself
function [probes snp_ids data]=load_cel_probes(CEL_file, CDFStruct, probesets)

% global CELStruct
% clear global CELStruct

if strcmpi(CEL_file(end-3:end),'.CEL')
    CELStruct=affyread(CEL_file);
elseif strcmpi(CEL_file(end-3:end),'.mat')
    load(CEL_file);
end

[dum loc]=ismember(probesets,{CDFStruct.ProbeSets.Name});
if ~all(loc)
    error('Not all probesets found');
end
probes=cell(1,56);
for i=1:7
    probes{(i-1)*8+1}=['Q' num2str(i) '_A_Sense_PM'];
    probes{(i-1)*8+2}=['Q' num2str(i) '_B_Sense_PM'];
    probes{(i-1)*8+3}=['Q' num2str(i) '_A_Sense_MM'];
    probes{(i-1)*8+4}=['Q' num2str(i) '_B_Sense_MM'];
    probes{(i-1)*8+5}=['Q' num2str(i) '_A_Antisense_PM'];
    probes{(i-1)*8+6}=['Q' num2str(i) '_B_Antisense_PM'];
    probes{(i-1)*8+7}=['Q' num2str(i) '_A_Antisense_MM'];
    probes{(i-1)*8+8}=['Q' num2str(i) '_B_Antisense_MM'];
end

data=nan(length(probesets),56,'single');

snp_ids = {};
for i=1:length(probesets)
    psvals = probesetvalues_no_bg(CELStruct,CDFStruct,loc(i));
    probeStruct = snpquartets(psvals,CDFStruct);
    snp_ids{i} = CDFStruct.ProbeSets(loc(i)).Name;
    
    for j=1:length(probeStruct.Quartet)
        if ~isempty(probeStruct.Quartet(j).A_Sense_PM)
            data(i,(j-1)*8+1)=probeStruct.Quartet(j).A_Sense_PM;
        end
        if ~isempty(probeStruct.Quartet(j).B_Sense_PM)
        data(i,(j-1)*8+2)=probeStruct.Quartet(j).B_Sense_PM;
        end
        if ~isempty(probeStruct.Quartet(j).A_Sense_MM)
        data(i,(j-1)*8+3)=probeStruct.Quartet(j).A_Sense_MM;
        end
        if ~isempty(probeStruct.Quartet(j).B_Sense_MM)
        data(i,(j-1)*8+4)=probeStruct.Quartet(j).B_Sense_MM;
        end
        if ~isempty(probeStruct.Quartet(j).A_Antisense_PM)
        data(i,(j-1)*8+5)=probeStruct.Quartet(j).A_Antisense_PM;
        end
        if ~isempty(probeStruct.Quartet(j).B_Antisense_PM)
        data(i,(j-1)*8+6)=probeStruct.Quartet(j).B_Antisense_PM;
        end
        if ~isempty(probeStruct.Quartet(j).A_Antisense_MM)
        data(i,(j-1)*8+7)=probeStruct.Quartet(j).A_Antisense_MM;
        end
        if ~isempty(probeStruct.Quartet(j).B_Antisense_MM)
        data(i,(j-1)*8+8)=probeStruct.Quartet(j).B_Antisense_MM;
        end
    end

    
    
    if mod(i,5000)==0
        disp(i);
    end
end

snp_ids = snp_ids';

