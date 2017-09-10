% Read the example files with 10 different SNPs from the CHIAMO algorithm
% homepage - these are the WCCC Welcome Trust samples.
function [A_copy B_copy] = ReadChiamoCohorts(chiamo_dir)

for c=1:9
    cohort{c} = loadCellFile(fullfile(chiamo_dir, ['Cohort' num2str(c) '.txt']));
    
    cmat{c} = char(cohort{c}(2:end,6:end));
    cmat{c} = vec_into_mat(str2num(cmat{c}),10);
    
    A_copy{c} = cmat{c}(:,1:2:end);     B_copy{c} = cmat{c}(:,2:2:end);
end



