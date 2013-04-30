% This function compares how well are the abberations based on two chips
% (e.g. xba and hind) for the same sample
%function [NumCopyChanges CopyNumCountsMat] = CompareTwoChipsAbberations(DispStruct, chroms)
function TwoDispStructsCopyChangesScore = CompareTwoDispStructsAbberations(DispStruct1, DispStruct2, chroms)

TwoDispStructsCopyChangesScore = 0;

TwoChipsTotalLength = 0; % how many SNPs do we have totally (used for normalizing the score)

% loop over chromosomes
for t=1:length(chroms)
    i=chroms(t);

    num_snps1 = length(DispStruct1.Chrom{i}.Locs);
    num_snps2 = length(DispStruct2.Chrom{i}.Locs);

    % combine both chips locations
    [joint_locations I] = sort([DispStruct1.Chrom{i}.Locs' DispStruct2.Chrom{i}.Locs']);
    I_inv = inv_perm(I);
    
    TotalCopyVec1 = zeros(1, length(joint_locations));
    TotalCopyVec2 = zeros(1, length(joint_locations));
        
    TotalCopyVec1(I_inv(1:num_snps1)) = DispStruct1.Chrom{i}.Data(:,1) + DispStruct1.Chrom{i}.Data(:,2);
    TotalCopyVec2(I_inv(num_snps1+1:end)) = DispStruct2.Chrom{i}.Data(:,1) + DispStruct2.Chrom{i}.Data(:,2);

    % perform interpolation
    TotalCopyVec1(I_inv(num_snps1+1:end)) = interp1(DispStruct1.Chrom{i}.Locs, TotalCopyVec1(I_inv(1:num_snps1)), DispStruct2.Chrom{i}.Locs);
    TotalCopyVec2(I_inv(1:num_snps1)) = interp1(DispStruct2.Chrom{i}.Locs, TotalCopyVec2(I_inv(num_snps1+1:end)), DispStruct1.Chrom{i}.Locs);


    % deal with NaNs at the boundaries:
    nan_inds1 = find(isnan(TotalCopyVec1));
    if(~isempty(nan_inds1))
        if(nan_inds1(1) == 1) % need to fill the start
            max_nan_inds = find(nan_inds1 == (1:length(nan_inds1)), 1, 'last' );
            TotalCopyVec1(1:max_nan_inds) = TotalCopyVec1(max_nan_inds+1);
        end
        if(nan_inds1(end) == length(TotalCopyVec1))
            max_nan_inds = find(nan_inds1 == ([1:length(nan_inds1)] + length(TotalCopyVec1) - length(nan_inds1) ), 1, 'first' );
            TotalCopyVec1(nan_inds1(max_nan_inds):end) = TotalCopyVec1(nan_inds1(max_nan_inds)-1);
        end
    end

    nan_inds2 = find(isnan(TotalCopyVec2));
    if(~isempty(nan_inds2))
        if(nan_inds2(1) == 1) % need to fill the start
            max_nan_inds = find(nan_inds2 == (1:length(nan_inds2)), 1, 'last' );
            TotalCopyVec2(1:max_nan_inds) = TotalCopyVec2(max_nan_inds+1);
        end
        if(nan_inds2(end) == length(TotalCopyVec2))
            max_nan_inds = find(nan_inds2 == ([1:length(nan_inds2)] + length(TotalCopyVec2) - length(nan_inds2) ), 1, 'first' );
            TotalCopyVec2(nan_inds2(max_nan_inds):end) = TotalCopyVec2(nan_inds2(max_nan_inds)-1);
        end
    end

    % Now compare the two vectors - what scoring function should we use?
    TwoDispStructsCopyChangesScore = TwoDispStructsCopyChangesScore + sum( (TotalCopyVec1 - TotalCopyVec2) .^ 2 );
    TwoChipsTotalLength = TwoChipsTotalLength + length(TotalCopyVec1);

end

TwoDispStructsCopyChangesScore = TwoDispStructsCopyChangesScore ./ TwoChipsTotalLength; % Normalization

