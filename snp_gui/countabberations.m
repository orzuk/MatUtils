% Count abberations for one sample
%function [NumCopyChanges CopyNumCountsMat] = CountAbberations(DispStruct, chroms)
function [NumCopyChanges CopyNumCountsMat] = CountAbberations(DispStruct, chroms)
    
NumCopyChanges = 0;
CopyNumCountsMat = zeros(1, 5);

% loop over chromosomes
for t=1:length(chroms)
    i=chroms(t);
    TotalCopyVec = DispStruct.Chrom{i}.Data(:,1) + DispStruct.Chrom{i}.Data(:,2);
    if(i ~= 23) % avoid chromosome X in copy # counts
        for j=1:5
            CopyNumCountsMat(j) = CopyNumCountsMat(j) + sum(TotalCopyVec == j-1);
        end
        NumCopyChanges = NumCopyChanges + sum(diff(TotalCopyVec) ~= 0);
    end
end