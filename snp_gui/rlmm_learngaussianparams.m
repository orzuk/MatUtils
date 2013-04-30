% Written by Or Zuk 7/2007
%
% Here we implement the RLMM package of Rabbee&Speed (Bioinformatics 2006),
% which was originally written in R, in Matlab.
% This file is part of the genotypng algorithm. It calculates the moments
% of the Multi-Dimensional Gaussian we shall use later in order to classify
% the raw SNP intensities into genotypes.
% We have some deviation from their algorithm, by taking a bayesian version
% ('BRLMM' by Affymetrix)
%
% The input to the function is a training set: matrices of A&B intensities
% and the corresponding correct genotypes.
%
% The input:
%
% CopyMatA - Intensities of A allele
% CopyMatB - Intensities of B allele
% GenotypesMat- Genotypes
% MinSparse - The minimum number such that above it the SNP is not considered sparase e.g. 5)
% PsuedoCount - The number of artificial 'samples' hypothetically used  to get the overall mean (e.g. 40)
%
% The output:
%
% RLMM Struct
% MuMats - matrix of Gaussian means
% Sigma_mats - matrix of Gaussian means
% SparseSNPs - For which SNPs we've used Regression since they were too sparse
function RLMM = RLMM_LearnGaussianParams(CopyMatA, CopyMatB, GenotypesMat, strand, chr, gender, MinSparse, PsuedoCount)

AssignAllGlobalConstants();

NumSNPs = size(CopyMatA,1)
NumSamples = size(CopyMatA,2)

% swap according to strand
minus_strand = strmatch('-', strand);
[CopyMatA(minus_strand,:) CopyMatB(minus_strand,:)] = swap(CopyMatA(minus_strand,:), CopyMatB(minus_strand,:));

% New: add to the sparse matrix also the SNPs of MALEs in the X chromosome - we learn the X's parameters only from the women
X_chr_inds = find(chr == 23); 
Male_inds = strmatch('M', gender);
CopyMatA(X_chr_inds, Male_inds) = NaN;
CopyMatB(X_chr_inds, Male_inds) = NaN;

% First find 'bad' SNPs - i.e. ones with 0 or 1 samples in one of the three groups
nan_mat = sparse(isnan(CopyMatA) | isnan(CopyMatB));

size(sum(GenotypesMat == AA,2))
size(sum(sparse(nan_mat & (GenotypesMat == AA))))




AA_Counts =  sum(GenotypesMat == AA,2) - sum(sparse(nan_mat & (GenotypesMat == AA)),2);
AB_Counts =  sum(GenotypesMat == AB,2) - sum(sparse(nan_mat & (GenotypesMat == AB)),2);
BB_Counts =  sum(GenotypesMat == BB,2) - sum(sparse(nan_mat & (GenotypesMat == BB)),2);

AA_CountMat = zeros(NumSNPs, NumSamples, 'uint8'); AA_CountMat(find(GenotypesMat == AA)) = 1;
AB_CountMat = zeros(NumSNPs, NumSamples, 'uint8'); AB_CountMat(find(GenotypesMat == AB)) = 1;
BB_CountMat = zeros(NumSNPs, NumSamples, 'uint8'); BB_CountMat(find(GenotypesMat == BB)) = 1;


% The Sigma matrix is as follows: 1 is A variance, 2 is B variance and 3 is the covariance
MAD = 1;
if(MAD)
    MuMats.AA = zeros(NumSNPs,2); MuMats.AB = zeros(NumSNPs,2); MuMats.BB = zeros(NumSNPs,2);
    SigmaMats.AA = zeros(NumSNPs,3); SigmaMats.AB = zeros(NumSNPs,3); SigmaMats.BB = zeros(NumSNPs,3);

    for i=1:NumSNPs
        AA_Inds = find(AA_CountMat(i,:));
        MuMats.AA(i,1) = nanmedian(CopyMatA(i,AA_Inds));
        MuMats.AA(i,2) = nanmedian(CopyMatB(i,AA_Inds));
        AB_Inds = find(AB_CountMat(i,:));
        MuMats.AB(i,1) = nanmedian(CopyMatA(i,AB_Inds));
        MuMats.AB(i,2) = nanmedian(CopyMatB(i,AB_Inds));
        BB_Inds = find(BB_CountMat(i,:));
        MuMats.BB(i,1) = nanmedian(CopyMatA(i,BB_Inds));
        MuMats.BB(i,2) = nanmedian(CopyMatB(i,BB_Inds));

        if(mod(i, 1000) == 0)
            mad_i_is = i
        end
    end

    % This consumes lots of memory - to remove ???
    %     DevMat.AA{1} = CopyMatA - repmat(MuMats.AA(:,1), 1, NumSamples);
    %     DevMat.AA{2} = CopyMatB - repmat(MuMats.AA(:,2), 1, NumSamples);
    %     DevMat.AB{1} = CopyMatA - repmat(MuMats.AB(:,1), 1, NumSamples);
    %     DevMat.AB{2} = CopyMatB - repmat(MuMats.AB(:,2), 1, NumSamples);
    %     DevMat.BB{1} = CopyMatA - repmat(MuMats.BB(:,1), 1, NumSamples);
    %     DevMat.BB{2} = CopyMatB - repmat(MuMats.BB(:,2), 1, NumSamples);

    for i=1:NumSNPs
        % This are one for every sample - may still contain NaNs
        DevMat.AA{1} = CopyMatA(i,:) - repmat(MuMats.AA(i,1), 1, NumSamples);
        DevMat.AA{2} = CopyMatB(i,:) - repmat(MuMats.AA(i,2), 1, NumSamples);
        DevMat.AB{1} = CopyMatA(i,:) - repmat(MuMats.AB(i,1), 1, NumSamples);
        DevMat.AB{2} = CopyMatB(i,:) - repmat(MuMats.AB(i,2), 1, NumSamples);
        DevMat.BB{1} = CopyMatA(i,:) - repmat(MuMats.BB(i,1), 1, NumSamples);
        DevMat.BB{2} = CopyMatB(i,:) - repmat(MuMats.BB(i,2), 1, NumSamples);


        AA_Inds = find(AA_CountMat(i,:));
        SigmaMats.AA(i,1) = MAD_CONST_SQR * nanmedian((DevMat.AA{1}(AA_Inds).^2));
        SigmaMats.AA(i,2) = MAD_CONST_SQR * nanmedian((DevMat.AA{2}(AA_Inds).^2));
        SigmaMats.AA(i,3) = MAD_CONST_SQR * nanmedian(DevMat.AA{1}(AA_Inds) .* DevMat.AA{2}(AA_Inds));   % Is this the same constant?
        AB_Inds = find(AB_CountMat(i,:));
        SigmaMats.AB(i,1) = MAD_CONST_SQR * nanmedian((DevMat.AB{1}(AB_Inds).^2));
        SigmaMats.AB(i,2) = MAD_CONST_SQR * nanmedian((DevMat.AB{2}(AB_Inds).^2));
        SigmaMats.AB(i,3) = MAD_CONST_SQR * nanmedian(DevMat.AB{1}(AB_Inds) .* DevMat.AB{2}(AB_Inds));   % Is this the same constant?
        BB_Inds = find(BB_CountMat(i,:));
        SigmaMats.BB(i,1) = MAD_CONST_SQR * nanmedian((DevMat.BB{1}(BB_Inds).^2));
        SigmaMats.BB(i,2) = MAD_CONST_SQR * nanmedian((DevMat.BB{2}(BB_Inds).^2));
        SigmaMats.BB(i,3) = MAD_CONST_SQR * nanmedian(DevMat.BB{1}(BB_Inds) .* DevMat.BB{2}(BB_Inds));   % Is this the same constant?

        % New correction: Use rotation, calculate variance and rotate back
        C = [SigmaMats.AA(i,1) SigmaMats.AA(i,3); SigmaMats.AA(i,3) SigmaMats.AA(i,2)];
        if(isempty(find(isnan(C), 1)))
            %%            if(~isposdef(C))
            [U,V] = eig(C);
            rotated_data = [DevMat.AA{1}(AA_Inds)' DevMat.AA{2}(AA_Inds)'] * U;
            W(1,1) = MAD_CONST_SQR .* mad(rotated_data(:,1),1).^2;
            W(2,2) = MAD_CONST_SQR .* mad(rotated_data(:,2),1).^2;
            W = min(MAX_VAR, max(W, MIN_VAR * eye(2)));
            C2 = U*W*inv(U); SigmaMats.AA(i,1) = C2(1,1);  SigmaMats.AA(i,2) = C2(2,2);  SigmaMats.AA(i,3) = C2(1,2);
            %%            end
        end
        C = [SigmaMats.AB(i,1) SigmaMats.AB(i,3); SigmaMats.AB(i,3) SigmaMats.AB(i,2)];
        if(isempty(find(isnan(C), 1)))
            %%           if(~isposdef(C))
            [U,V] = eig(C);
            rotated_data = [DevMat.AB{1}(AB_Inds)' DevMat.AB{2}(AB_Inds)'] * U;
            W(1,1) = MAD_CONST_SQR .* mad(rotated_data(:,1),1).^2;
            W(2,2) = MAD_CONST_SQR .* mad(rotated_data(:,2),1).^2;
            W = min(MAX_VAR, max(W, MIN_VAR * eye(2)));
            C2 = U*W*inv(U); SigmaMats.AB(i,1) = C2(1,1);  SigmaMats.AB(i,2) = C2(2,2);  SigmaMats.AB(i,3) = C2(1,2);
            %%            end
        end
        C = [SigmaMats.BB(i,1) SigmaMats.BB(i,3); SigmaMats.BB(i,3) SigmaMats.BB(i,2)];
        if(isempty(find(isnan(C), 1)))
            %%            if(~isposdef(C))
            [U,V] = eig(C);
            rotated_data = [DevMat.BB{1}(BB_Inds)' DevMat.BB{2}(BB_Inds)'] * U;
            W(1,1) = MAD_CONST_SQR .* mad(rotated_data(:,1),1).^2;
            W(2,2) = MAD_CONST_SQR .* mad(rotated_data(:,2),1).^2;
            W = min(MAX_VAR, max(W, MIN_VAR * eye(2)));
            C2 = U*W*inv(U); SigmaMats.BB(i,1) = C2(1,1);  SigmaMats.BB(i,2) = C2(2,2);  SigmaMats.BB(i,3) = C2(1,2);
            %%            end
        end

        if(mod(i, 1000) == 0)
            learn_snp_i_is = i
        end

    end


    % Correction to force matrices to be positive-definite
    %     SigmaMats.AA(:,3) = max( -sqrt(SigmaMats.AA(:,1).*SigmaMats.AA(:,2))+EPSILON, ...
    %         min( SigmaMats.AA(:,3), sqrt(SigmaMats.AA(:,1).*SigmaMats.AA(:,2))-EPSILON ) );
    %     SigmaMats.AB(:,3) = max( -sqrt(SigmaMats.AB(:,1).*SigmaMats.AB(:,2))+EPSILON, ...
    %         min( SigmaMats.AB(:,3), sqrt(SigmaMats.AB(:,1).*SigmaMats.AB(:,2))-EPSILON ) );
    %     SigmaMats.BB(:,3) = max( -sqrt(SigmaMats.BB(:,1).*SigmaMats.BB(:,2))+EPSILON, ...
    %         min( SigmaMats.BB(:,3), sqrt(SigmaMats.BB(:,1).*SigmaMats.BB(:,2))-EPSILON ) );


else  % Probably not working because of NaNs ...
    MuMats.AA(:,1) = sum(CopyMatA .* single(AA_CountMat),2) ./ AA_Counts;  % Mean A
    MuMats.AA(:,2) = sum(CopyMatB .* single(AA_CountMat),2) ./ AA_Counts;  % Mean B
    MuMats.AB(:,1) = sum(CopyMatA .* single(AB_CountMat),2) ./ AB_Counts;
    MuMats.AB(:,2) = sum(CopyMatB .* single(AB_CountMat),2) ./ AB_Counts;
    MuMats.BB(:,1) = sum(CopyMatA .* single(BB_CountMat),2) ./ BB_Counts;
    MuMats.BB(:,2) = sum(CopyMatB .* single(BB_CountMat),2) ./ BB_Counts;

    SigmaMats.AA(:,1) = sum((CopyMatA.^2) .* single(AA_CountMat), 2) ./ AA_Counts - MuMats.AA(:,1).^2; % Var A
    SigmaMats.AA(:,2) = sum((CopyMatB.^2) .* single(AA_CountMat), 2) ./ AA_Counts - MuMats.AA(:,2).^2; % Var B
    SigmaMats.AA(:,3) = sum(CopyMatA .* CopyMatB .* single(AA_CountMat), 2) ./ AA_Counts - ...
        MuMats.AA(:,1) .* MuMats.AA(:,2);   % Cov(A,B)
    SigmaMats.AB(:,1) = sum((CopyMatA.^2) .* single(AB_CountMat), 2) ./ AB_Counts - MuMats.AB(:,1).^2; % Var A
    SigmaMats.AB(:,2) = sum((CopyMatB.^2) .* single(AB_CountMat), 2) ./ AB_Counts - MuMats.AB(:,2).^2; % Var B
    SigmaMats.AB(:,3) = sum(CopyMatA .* CopyMatB .* single(AB_CountMat), 2) ./ AB_Counts - ...
        MuMats.AB(:,1) .* MuMats.AB(:,2);   % Cov(A,B)
    SigmaMats.BB(:,1) = sum((CopyMatA.^2) .* single(BB_CountMat), 2) ./ BB_Counts - MuMats.BB(:,1).^2; % Var A
    SigmaMats.BB(:,2) = sum((CopyMatB.^2) .* single(BB_CountMat), 2) ./ BB_Counts - MuMats.BB(:,2).^2; % Var B
    SigmaMats.BB(:,3) = sum(CopyMatA .* CopyMatB .* single(BB_CountMat), 2) ./ BB_Counts - ...
        MuMats.BB(:,1) .* MuMats.BB(:,2);   % Cov(A,B)
end

MinSparse = max(MinSparse, 2); % We do not allow one sample to work on
SparseSNPsAA = find(AA_Counts < MinSparse);
SparseSNPsAB = find(AB_Counts < MinSparse);
SparseSNPsBB = find(BB_Counts < MinSparse);
AllSparseSNPs = intersect( intersect(SparseSNPsAA, SparseSNPsAB),  SparseSNPsBB ); % something wrong - no calls for these SNPs
TwoSparseSNPs = find(median([AA_Counts AB_Counts BB_Counts],2) < MinSparse); % In these SNPs there's only one full allele
SparseSNPs = find(min( min(AA_Counts, AB_Counts), BB_Counts) < MinSparse);  % In these SNPs there's one sparse allele
OneSparseSNPs = setdiff(SparseSNPs, TwoSparseSNPs); TwoSparseSNPs = setdiff(TwoSparseSNPs, AllSparseSNPs);
GoodSNPs = setdiff(1:NumSNPs, SparseSNPs);

OneSparseCounts = [AA_Counts(OneSparseSNPs) AB_Counts(OneSparseSNPs) BB_Counts(OneSparseSNPs)];
[OneMinCounts OneMinInds] = min(OneSparseCounts,[],2); % Here we take the only sparse index
TwoSparseCounts = [AA_Counts(TwoSparseSNPs) AB_Counts(TwoSparseSNPs) BB_Counts(TwoSparseSNPs)];
[TwoMaxCounts TwoMaxInds] = max(TwoSparseCounts,[],2); % Here we take the only non-sparse index



MuMuVec = [ mean(MuMats.AA(GoodSNPs,:) )  mean(MuMats.AB(GoodSNPs,:) )  mean(MuMats.BB(GoodSNPs,:) ) ];
MuSigmaMat = cov(  [MuMats.AA(GoodSNPs,:) MuMats.AB(GoodSNPs,:) MuMats.BB(GoodSNPs,:)] ); % corr
% Get the partial matrices for the conditional distributions
PartialMuMuVec{AA}.m1 = MuMuVec(1:2); PartialMuMuVec{AA}.m2 = MuMuVec(3:6);
PartialMuSigmaMats{AA}.s12 = MuSigmaMat(1:2,3:6);
PartialMuSigmaMats{AA}.s22_Inv = inv(MuSigmaMat(3:6,3:6));
PartialMuSigmaMats{AA}.two_s12 = MuSigmaMat(3:6,1:2);
PartialMuSigmaMats{AA}.two_s22_Inv = inv(MuSigmaMat(1:2,1:2));

PartialMuMuVec{AB}.m1 = MuMuVec(3:4); PartialMuMuVec{AB}.m2 = MuMuVec([1:2,5:6]);
PartialMuSigmaMats{AB}.s12 = MuSigmaMat(3:4,[1:2,5:6]);
PartialMuSigmaMats{AB}.s22_Inv = inv(MuSigmaMat([1:2,5:6],[1:2,5:6]));
PartialMuSigmaMats{AB}.two_s12 = MuSigmaMat([1:2,5:6],3:4);
PartialMuSigmaMats{AB}.two_s22_Inv = inv(MuSigmaMat(3:4,3:4));

PartialMuMuVec{BB}.m1 = MuMuVec(5:6); PartialMuMuVec{BB}.m2 = MuMuVec(1:4);
PartialMuSigmaMats{BB}.s12 = MuSigmaMat(5:6,1:4);
PartialMuSigmaMats{BB}.s22_Inv = inv(MuSigmaMat(1:4,1:4));
PartialMuSigmaMats{BB}.two_s12 = MuSigmaMat(1:4,5:6);
PartialMuSigmaMats{BB}.two_s22_Inv = inv(MuSigmaMat(5:6,5:6));


SigmaMuVec = [mean(SigmaMats.AA(GoodSNPs,:)) mean(SigmaMats.AB(GoodSNPs,:)) mean(SigmaMats.BB(GoodSNPs,:))];
SigmaSigmaMat = cov( [ SigmaMats.AA(GoodSNPs,:) SigmaMats.AB(GoodSNPs,:) SigmaMats.BB(GoodSNPs,:)]); % corr

PartialSigmaMuVec{AA}.m1 = SigmaMuVec(1:3); PartialSigmaMuVec{AA}.m2 = SigmaMuVec(4:9);
PartialSigmaSigmaMats{AA}.s12 = SigmaSigmaMat(1:3,4:9);
PartialSigmaSigmaMats{AA}.s22_Inv = inv(SigmaSigmaMat(4:9,4:9));
PartialSigmaSigmaMats{AA}.two_s12 = SigmaSigmaMat(4:9,1:3);
PartialSigmaSigmaMats{AA}.two_s22_Inv = inv(SigmaSigmaMat(1:3,1:3));


PartialSigmaMuVec{AB}.m1 = SigmaMuVec(4:6); PartialSigmaMuVec{AB}.m2 = SigmaMuVec([1:3,7:9]);
PartialSigmaSigmaMats{AB}.s12 = SigmaSigmaMat(4:6,[1:3,7:9]);
PartialSigmaSigmaMats{AB}.s22_Inv = inv(SigmaSigmaMat([1:3,7:9],[1:3,7:9]));
PartialSigmaSigmaMats{AB}.two_s12 = SigmaSigmaMat([1:3,7:9],4:6);
PartialSigmaSigmaMats{AB}.two_s22_Inv = inv(SigmaSigmaMat(4:6,4:6));


PartialSigmaMuVec{BB}.m1 = SigmaMuVec(7:9); PartialSigmaMuVec{BB}.m2 = SigmaMuVec(1:6);
PartialSigmaSigmaMats{BB}.s12 = SigmaSigmaMat(7:9,1:6);
PartialSigmaSigmaMats{BB}.s22_Inv = inv(SigmaSigmaMat(1:6,1:6));
PartialSigmaSigmaMats{BB}.two_s12 = SigmaSigmaMat(1:6,7:9);
PartialSigmaSigmaMats{BB}.two_s22_Inv = inv(SigmaSigmaMat(7:9,7:9));


%% SavedMuMats = MuMats; SavedSigmaMats = SigmaMats;
% Find 'revesed' SNPs
ReverseSNPsAA = find(MuMats.AA(:,1) < MuMats.AA(:,2));
ReverseSNPsBB = find(MuMats.BB(:,1) > MuMats.BB(:,2));
reverse_AA_non_sparse = setdiff(ReverseSNPsAA, SparseSNPs);
ReverseSNPs = intersect(ReverseSNPsAA, ReverseSNPsBB);

% swap the SNPs with reverse indices
if(~isempty(ReverseSNPs))
    [CopyMatA(ReverseSNPs,:) CopyMatB(ReverseSNPs,:)] = swap(CopyMatA(ReverseSNPs,:), CopyMatB(ReverseSNPs,:));
    [MuMats.AA(ReverseSNPs,:) MuMats.BB(ReverseSNPs,:)] = swap(MuMats.AA(ReverseSNPs,:), MuMats.BB(ReverseSNPs,:));
    [SigmaMats.AA(ReverseSNPs,:) SigmaMats.BB(ReverseSNPs,:)] = swap(SigmaMats.AA(ReverseSNPs,:), SigmaMats.BB(ReverseSNPs,:));
end

ReverseSNPsAA = setdiff(ReverseSNPsAA, ReverseSNPs); 
ReverseSNPsBB = setdiff(ReverseSNPsBB, ReverseSNPs); 
if(~isempty(ReverseSNPsAA))
    MuMats.AA(ReverseSNPsAA,:) = repmat(MuMuVec(1:2), length(ReverseSNPsAA), 1);
    SigmaMats.AA(ReverseSNPsAA,:) = repmat(SigmaMuVec(1:3), length(ReverseSNPsAA), 1);
end
if(~isempty(ReverseSNPsBB))
    MuMats.BB(ReverseSNPsBB,:) = repmat(MuMuVec(5:6), length(ReverseSNPsBB), 1);
    SigmaMats.BB(ReverseSNPsBB,:) = repmat(SigmaMuVec(7:9), length(ReverseSNPsBB), 1);
end


ctr=1;
% First deal with the 'ALL SPARSE' SNPs (everything NN or something) - stays the same with BRLMM
if(~isempty(AllSparseSNPs))
    NumSparse = length(AllSparseSNPs);

    MuMats.AA(AllSparseSNPs,:) = repmat(MuMuVec(1:2), NumSparse, 1);
    MuMats.AB(AllSparseSNPs,:) = repmat(MuMuVec(3:4), NumSparse, 1);
    MuMats.BB(AllSparseSNPs,:) = repmat(MuMuVec(5:6), NumSparse, 1);

    SigmaMats.AA(AllSparseSNPs,:) = repmat(SigmaMuVec(1:3), NumSparse, 1);
    SigmaMats.AB(AllSparseSNPs,:) = repmat(SigmaMuVec(4:6), NumSparse, 1);
    SigmaMats.BB(AllSparseSNPs,:) = repmat(SigmaMuVec(7:9), NumSparse, 1);

end
% Here deal with the 'TWO SPARSE' SNPs (only one allele is present in high numbers)
if(~isempty(TwoSparseSNPs))
    % Currently don't do anything
    NumSparse = length(TwoSparseSNPs);
    % First calculate the a vectors
    MuA = [MuMats.AA(TwoSparseSNPs,:) MuMats.AB(TwoSparseSNPs,:) MuMats.BB(TwoSparseSNPs,:)]';
    J = [2*TwoMaxInds'-1 2*TwoMaxInds']; I = [1:NumSparse 1:NumSparse];
    IndsMat = sparse(I,J,ones(1,2*NumSparse),NumSparse,6); MuA = reshape(MuA(logical(IndsMat')), 2, NumSparse)';

    SigmaA = [SigmaMats.AA(TwoSparseSNPs,:) SigmaMats.AB(TwoSparseSNPs,:) SigmaMats.BB(TwoSparseSNPs,:)]';
    J = [3*TwoMaxInds'-2 3*TwoMaxInds'-1 3*TwoMaxInds']; I = [1:NumSparse 1:NumSparse 1:NumSparse];
    IndsMat = sparse(I,J,ones(1,3*NumSparse),NumSparse,9); SigmaA = reshape(SigmaA(logical(IndsMat')), 3, NumSparse)';
    % Note: The loop should be replaced in the future by a vectorized version !!!
    for i=TwoSparseSNPs' % Here we do not have enough data to learn from and need to do something else ..
        MuBar = PartialMuMuVec{TwoMaxInds(ctr)}.m2' + PartialMuSigmaMats{OneMinInds(ctr)}.two_s12 * PartialMuSigmaMats{OneMinInds(ctr)}.two_s22_Inv * ...
            (MuA(ctr,:) - PartialMuMuVec{TwoMaxInds(ctr)}.m1)';
        SigmaBar = PartialSigmaMuVec{TwoMaxInds(ctr)}.m2' + PartialSigmaSigmaMats{OneMinInds(ctr)}.two_s12 * PartialSigmaSigmaMats{OneMinInds(ctr)}.two_s22_Inv * ...
            (SigmaA(ctr,:) - PartialSigmaMuVec{TwoMaxInds(ctr)}.m1)';
        if(TwoMaxInds(ctr) == AA)
            if(AB_Counts(i) < 1)
                MuMats.AB(i,:) = MuBar(1:2)';
            else
                MuMats.AB(i,:) = (MuBar(1:2)' * PsuedoCount + MuMats.AB(i,:) * AB_Counts(i)) ./ (PsuedoCount + AB_Counts(i));
            end
            if(BB_Counts(i) < 1)
                MuMats.BB(i,:) = MuBar(3:4)';
            else
                MuMats.BB(i,:) = (MuBar(3:4)' * PsuedoCount + MuMats.BB(i,:) * BB_Counts(i)) ./ (PsuedoCount + BB_Counts(i));
            end

        end
        if(TwoMaxInds(ctr) == AB)
            if(AA_Counts(i) < 1)
                MuMats.AA(i,:) = MuBar(1:2)';
            else
                MuMats.AA(i,:) = (MuBar(1:2)' * PsuedoCount + MuMats.AA(i,:) * AA_Counts(i)) ./ (PsuedoCount + AA_Counts(i));
            end
            if(BB_Counts(i) < 1)
                MuMats.BB(i,:) = MuBar(3:4)';
            else
                MuMats.BB(i,:) = (MuBar(3:4)' * PsuedoCount + MuMats.BB(i,:) * BB_Counts(i)) ./ (PsuedoCount + BB_Counts(i));
            end
        end
        if(TwoMaxInds(ctr) == BB)
            if(AA_Counts(i) < 1)
                MuMats.AA(i,:) = MuBar(1:2)';
            else
                MuMats.AA(i,:) = (MuBar(1:2)' * PsuedoCount + MuMats.AA(i,:) * AA_Counts(i)) ./ (PsuedoCount + AA_Counts(i));
            end
            if(AB_Counts(i) < 1)
                MuMats.AB(i,:) = MuBar(3:4)';
            else
                MuMats.AB(i,:) = (MuBar(3:4)' * PsuedoCount + MuMats.AB(i,:) * AB_Counts(i)) ./ (PsuedoCount + AB_Counts(i));
            end
        end


        %%% New Bayesian Updating rule for the Sigmas (BRLMM)
        if(AA_Counts(i) < 2)
            SigmaMats.AA(i,:) = SigmaMuVec(1:3);
        else
            SigmaMats.AA(i,:) = (SigmaMuVec(1:3) * PsuedoCount + SigmaMats.AA(i,:) * (AA_Counts(i)-1)) ./  ...
                (PsuedoCount + AA_Counts(i)-1);
        end
        if(AB_Counts(i) < 2)
            SigmaMats.AB(i,:) = SigmaMuVec(4:6);
        else
            SigmaMats.AB(i,:) = (SigmaMuVec(4:6) * PsuedoCount + SigmaMats.AB(i,:) * (AB_Counts(i)-1)) ./  ...
                (PsuedoCount + AB_Counts(i)-1);
        end
        if(BB_Counts(i) < 2)
            SigmaMats.BB(i,:) = SigmaMuVec(7:9);
        else
            SigmaMats.BB(i,:) = (SigmaMuVec(7:9) * PsuedoCount + SigmaMats.BB(i,:) * (BB_Counts(i)-1)) ./  ...
                (PsuedoCount + BB_Counts(i)-1);
        end

        ctr=ctr+1;
    end
end
% Here deal with the 'ONE SPARSE' SNPs
ctr=1;
if(~isempty(OneSparseSNPs))
    NumSparse = length(OneSparseSNPs);
    % First calculate the a vectors
    MuA = [MuMats.AA(OneSparseSNPs,:) MuMats.AB(OneSparseSNPs,:) MuMats.BB(OneSparseSNPs,:)]';
    J = [2*OneMinInds'-1 2*OneMinInds']; I = [1:NumSparse 1:NumSparse];
    IndsMat = sparse(I,J,ones(1,2*NumSparse),NumSparse,6); MuA = reshape(MuA(logical(1-IndsMat')), 4, NumSparse)';

    SigmaA = [SigmaMats.AA(OneSparseSNPs,:) SigmaMats.AB(OneSparseSNPs,:) SigmaMats.BB(OneSparseSNPs,:)]';
    J = [3*OneMinInds'-2 3*OneMinInds'-1 3*OneMinInds']; I = [1:NumSparse 1:NumSparse 1:NumSparse];
    IndsMat = sparse(I,J,ones(1,3*NumSparse),NumSparse,9); SigmaA = reshape(SigmaA(logical(1-IndsMat')), 6, NumSparse)';
    % Note: The loop should be replaced in the future by a vectorized version !!!
    for i=OneSparseSNPs' % Here we do not have enough data to learn from and need to do something else ..
        MuBar = PartialMuMuVec{OneMinInds(ctr)}.m1' + PartialMuSigmaMats{OneMinInds(ctr)}.s12 * PartialMuSigmaMats{OneMinInds(ctr)}.s22_Inv * ...
            (MuA(ctr,:) - PartialMuMuVec{OneMinInds(ctr)}.m2)';
        SigmaBar = PartialSigmaMuVec{OneMinInds(ctr)}.m1' + PartialSigmaSigmaMats{OneMinInds(ctr)}.s12 * PartialSigmaSigmaMats{OneMinInds(ctr)}.s22_Inv * ...
            (SigmaA(ctr,:) - PartialSigmaMuVec{OneMinInds(ctr)}.m2)';


        if(OneMinInds(ctr) == AA)
            if(AA_Counts(i) < 1)
                MuMats.AA(i,:) = MuBar';
            else
                MuMats.AA(i,:) = (MuBar' * PsuedoCount + MuMats.AA(i,:) * AA_Counts(i)) ./ (PsuedoCount + AA_Counts(i));
            end
        end
        if(OneMinInds(ctr) == AB)
            if(AB_Counts(i) < 1)
                MuMats.AB(i,:) = MuBar';
            else
                MuMats.AB(i,:) = (MuBar' * PsuedoCount + MuMats.AB(i,:) * AB_Counts(i)) ./ (PsuedoCount + AB_Counts(i));
            end
        end
        if(OneMinInds(ctr) == BB)
            if(BB_Counts(i) < 1)
                MuMats.BB(i,:) = MuBar';
            else
                MuMats.BB(i,:) = (MuBar' * PsuedoCount + MuMats.BB(i,:) * BB_Counts(i)) ./ (PsuedoCount + BB_Counts(i));
            end
        end

        %%% New Bayesian Updating rule for the Sigmas (BRLMM)
        if(AA_Counts(i) < 2)
            SigmaMats.AA(i,:) = SigmaMuVec(1:3);
        else
            SigmaMats.AA(i,:) = (SigmaMuVec(1:3) * PsuedoCount + SigmaMats.AA(i,:) * (AA_Counts(i)-1)) ./  ...
                (PsuedoCount + AA_Counts(i)-1);
        end
        if(AB_Counts(i) < 2)
            SigmaMats.AB(i,:) = SigmaMuVec(4:6);
        else
            SigmaMats.AB(i,:) = (SigmaMuVec(4:6) * PsuedoCount + SigmaMats.AB(i,:) * (AB_Counts(i)-1)) ./  ...
                (PsuedoCount + AB_Counts(i)-1);
        end
        if(BB_Counts(i) < 2)
            SigmaMats.BB(i,:) = SigmaMuVec(7:9);
        else
            SigmaMats.BB(i,:) = (SigmaMuVec(7:9) * PsuedoCount + SigmaMats.BB(i,:) * (BB_Counts(i)-1)) ./  ...
                (PsuedoCount + BB_Counts(i)-1);
        end

        ctr=ctr+1;
    end
end



%%     % Plot for comparison between the 'sparse' and the 'dense' estimations - MU
%%     figure; hold on; plot(mat_into_vec(SavedMuMats.AA(SparseSNPs,:)), mat_into_vec(MuMats.AA(SparseSNPs,:)), '.');
%%     plot(mat_into_vec(SavedMuMats.AB(SparseSNPs,:)), mat_into_vec(MuMats.AB(SparseSNPs,:)), 'r.');
%%     plot(mat_into_vec(SavedMuMats.BB(SparseSNPs,:)), mat_into_vec(MuMats.BB(SparseSNPs,:)), 'g.');
%%     legend('AA', 'AB', 'BB'); title('Compare Dense to Sparse \mu estimation'); xlabel('Dense'); ylabel('Sparse');
%%
%%     % Plot for comparison between the 'sparse' and the 'dense' estimations - SIGMA
%%     figure; hold on; plot(mat_into_vec(SavedSi% Set the function's output
RLMM.MuMuVec = MuMuVec;
RLMM.MuSigmaMat = MuSigmaMat;
RLMM.SigmaMuVec = SigmaMuVec;
RLMM.SigmaSigmaMat = SigmaSigmaMat;
RLMM.MuMats = MuMats;
RLMM.SigmaMats = SigmaMats;
RLMM.SparseSNPs = SparseSNPs;
RLMM.GoodSNPs = GoodSNPs'; % Transpose to make orientations consistent.
RLMM.AllSparseSNPs = AllSparseSNPs;
RLMM.TwoSparseSNPs = TwoSparseSNPs;
RLMM.OneSparseSNPs = OneSparseSNPs;
RLMM.ReverseSNPsAA = ReverseSNPsAA; RLMM.ReverseSNPsBB = ReverseSNPsBB;
RLMM.ReverseSNPs = union(ReverseSNPsAA, ReverseSNPsBB); % Check which SNPs were flipped

%%     plot(mat_into_vec(SavedSigmaMats.AB(SparseSNPs,:)), mat_into_vec(SigmaMats.AB(SparseSNPs,:)), 'r.');
%%     plot(mat_into_vec(SavedSigmaMats.BB(SparseSNPs,:)), mat_into_vec(SigmaMats.BB(SparseSNPs,:)), 'g.');
%%     legend('AA', 'AB', 'BB'); title('Compare Dense to Sparse \Sigma estimation'); xlabel('Dense'); ylabel('Sparse')


%% MuMats = SavedMuMats; SigmaMats = SavedSigmaMats; % Get the dense ones back

