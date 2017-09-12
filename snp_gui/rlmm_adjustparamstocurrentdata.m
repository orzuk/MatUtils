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
function RLMM_Correction = RLMM_AdjustParamsToCurrentData(CopyMatA, CopyMatB, snp_ids, hapmap_population, ...
    chip_type, PsuedoCount)

AssignAllGlobalConstants();
NumSNPs = length(CopyMatA);
load(['../database/RLMM_' pop_str_vec{hapmap_population} '_' chip_type '.mat']);

for i=1:3
    INIT_S{i} = zeros(2);
    INIT_S{i}(1,1) = RLMM.SigmaMuVec((i-1)*3+1);
    INIT_S{i}(2,2) = RLMM.SigmaMuVec((i-1)*3+2);
    INIT_S{i}(1,2) = RLMM.SigmaMuVec((i-1)*3+3); INIT_S{i}(2,1) = INIT_S{i}(1,2);
end

[P,M,S, LogLike]=MixtureOfGaussiansMultiDimGivenInit([CopyMatA CopyMatB]',3,10, ...
     [0.375 0.375 0.25], vec_into_mat(RLMM.MuMuVec, 2)', INIT_S);

figure; plot(CopyMatA, CopyMatB, '.g');

labels_vec{1} ='xxx';
labels_vec{2} ='yyy';
legends_vec = {'AA', 'AB', 'BB'};

MixtureOfGaussiansDraw2dGaussians(vec_into_mat(RLMM.MuMuVec, 2)', INIT_S, labels_vec, legends_vec)
MixtureOfGaussiansDraw2dGaussians(M, S, labels_vec, legends_vec)

[S I J] = intersect(snp_ids, RLMM.snp_ids);
MuMats.AA = repmat(RLMM.MuMuVec(1:2), NumSNPs, 1);
MuMats.AB = repmat(RLMM.MuMuVec(3:4), NumSNPs, 1);
MuMats.BB = repmat(RLMM.MuMuVec(5:6), NumSNPs, 1);
SigmaMats.AA = repmat(RLMM.SigmaMuVec(1:3), NumSNPs, 1);
SigmaMats.AB = repmat(RLMM.SigmaMuVec(4:6), NumSNPs, 1);
SigmaMats.BB = repmat(RLMM.SigmaMuVec(7:9), NumSNPs, 1);
MuMats.AA(I,:) = RLMM.MuMats.AA(J,:);
MuMats.AB(I,:) = RLMM.MuMats.AB(J,:);
MuMats.BB(I,:) = RLMM.MuMats.BB(J,:);
SigmaMats.AA(I,:) = RLMM.SigmaMats.AA(J,:);
SigmaMats.AB(I,:) = RLMM.SigmaMats.AB(J,:);
SigmaMats.BB(I,:) = RLMM.SigmaMats.BB(J,:);

max_vec = max( max( max(SigmaMats.AA(:,1:2)'), max(SigmaMats.AB(:,1:2)')), max(SigmaMats.BB(:,1:2)'));

[LearnedGenotypesMat Mahalahonis] = RLMM_ClassifyUsingGaussianParams(CopyMatA, CopyMatB, MuMats, SigmaMats);



ADDITIVE_CORRECTION = 0; MULTIPLICATIVE_CORRECTION = 1;
correction_type = ADDITIVE_CORRECTION;
if(correction_type == ADDITIVE_CORRECTION)
    AdjustMuMats.AA = MuMats.AA + repmat( M(1,:)-RLMM.MuMuVec(1:2), NumSNPs, 1);
    AdjustMuMats.AB = MuMats.AB + repmat( M(2,:)-RLMM.MuMuVec(3:4), NumSNPs, 1);
    AdjustMuMats.BB = MuMats.BB + repmat( M(3,:)-RLMM.MuMuVec(5:6), NumSNPs, 1);
else % here MULTIPLICATIVE_CORRECTION
    AdjustMuMats.AA = MuMats.AA .* repmat( M(1,:)./RLMM.MuMuVec(1:2), NumSNPs, 1);
    AdjustMuMats.AB = MuMats.AB .* repmat( M(1,:)./RLMM.MuMuVec(1:2), NumSNPs, 1);
    AdjustMuMats.BB = MuMats.BB .* repmat( M(1,:)./RLMM.MuMuVec(1:2), NumSNPs, 1);
end

[AdjustLearnedGenotypesMat AdjustMahalahonis] = RLMM_ClassifyUsingGaussianParams(CopyMatA, CopyMatB, ...
    AdjustMuMats, SigmaMats);

PlotAlleleRatios(CopyMatA+CopyMatB, CopyMatA./CopyMatB, LearnedGenotypesMat, 'Hapmap RLMM Params');
PlotAlleleRatios(CopyMatA+CopyMatB, CopyMatA./CopyMatB, AdjustLearnedGenotypesMat, 'Adjusted RLMM Params');
T=load('F:\Or\LeukemiaDchip\HD78_9_n_xba.mat'); % load old -dchip Eytan14-pc
%% T=load('E:\public\SNP_GUI\SNP_GUI\data\LeukemiaDchip\HD78_9_n_xba.mat'); % load old -dchip my pc
[SSS III JJJ] = intersect(snp_ids, T.snp_id_xba);
strand = load_snp_ids_strand(T.snp_id_xba, chip_type);
% change genotype according to strand
minus_strand = strmatch('-', strand);

temp = T.genotype_vec_xba;
T.genotype_vec_xba(intersect(find(temp==AA),minus_strand)) = BB;
T.genotype_vec_xba(intersect(find(temp==BB),minus_strand)) = AA;

Errors = find(LearnedGenotypesMat(III) ~= T.genotype_vec_xba(JJJ) & (T.genotype_vec_xba(JJJ) ~= 4) );
AdjustErrors = find( (AdjustLearnedGenotypesMat(III) ~= T.genotype_vec_xba(JJJ)) & (T.genotype_vec_xba(JJJ) ~= 4) );
AdjustErrors2 = sum(AdjustLearnedGenotypesMat(III) ~= LearnedGenotypesMat(III))

Errors_AA = intersect(Errors, find(LearnedGenotypesMat(III) == AA));
Errors_AB = intersect(Errors, find(LearnedGenotypesMat(III) == AB));
Errors_BB = intersect(Errors, find(LearnedGenotypesMat(III) == BB));


Errors_True_AA = intersect(Errors, find(T.genotype_vec_xba(JJJ) == AA));
Errors_True_AB = intersect(Errors, find(T.genotype_vec_xba(JJJ) == AB));
Errors_True_BB = intersect(Errors, find(T.genotype_vec_xba(JJJ) == BB));


PlotAlleleRatios(CopyMatA(III(Errors))+CopyMatB(III(Errors)), CopyMatA(III(Errors))./CopyMatB(III(Errors)), ...
    LearnedGenotypesMat(III(Errors)), 'BAD SNPs Hapmap RLMM Params');
PlotAlleleRatios(CopyMatA(III(AdjustErrors))+CopyMatB(III(AdjustErrors)), CopyMatA(III(AdjustErrors))./CopyMatB(III(AdjustErrors)), ...
    AdjustLearnedGenotypesMat(III(AdjustErrors)), 'BAD SNPs Adjusted RLMM Params');

PlotAlleleRatios(CopyMatA(III(Errors))+CopyMatB(III(Errors)), CopyMatA(III(Errors))./CopyMatB(III(Errors)), ...
    T.genotype_vec_xba(JJJ(Errors)), 'AFFY BAD SNPs Hapmap RLMM Params');
PlotAlleleRatios(CopyMatA(III(AdjustErrors))+CopyMatB(III(AdjustErrors)), CopyMatA(III(AdjustErrors))./CopyMatB(III(AdjustErrors)), ...
    T.genotype_vec_xba(JJJ(AdjustErrors)), 'AFFY BAD SNPs Adjusted RLMM Params');

figure; hold on; plot(MuMats.AA(III(Errors),1), MuMats.AA(III(Errors),2), '.');
plot(MuMats.AB(III(Errors),1), MuMats.AB(III(Errors),2), 'g.');
plot(MuMats.BB(III(Errors),1), MuMats.BB(III(Errors),2), '.r');

figure; hold on; plot(MuMats.AA(:,1), MuMats.AA(:,2), '.');
plot(MuMats.AB(:,1), MuMats.AB(:,2), 'g.');
plot(MuMats.BB(:,1), MuMats.BB(:,2), '.r');



XXX = III(Errors);


[S_ERR I_ERR J_ERR] = intersect(snp_ids(XXX), RLMM.snp_ids);

ALL_DATA = load('E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\CEU\AllSamplesMat_xba.mat');

[S_ALL I_ALL J_ALL] = intersect(snp_ids(x), ALL_DATA.snp_ids);
figure; plot(ALL_DATA.NormalizedSNPsCopyMatA(J_ALL,:), ALL_DATA.NormalizedSNPsCopyMatB(J_ALL,:), '.');

x=III(Errors(5))
%% PlotAlleleRatios(CopyMatA(x)+CopyMatB(x), CopyMatA(x)./CopyMatB(x), LearnedGenotypesMat(x), 'BAD SNPs Hapmap RLMM Params');
figure;
SIG{1} = [SigmaMats.AA(x,1) SigmaMats.AA(x,3); SigmaMats.AA(x,3) SigmaMats.AA(x,2)]; 
SIG{2} = [SigmaMats.AB(x,1) SigmaMats.AB(x,3); SigmaMats.AB(x,3) SigmaMats.AB(x,2)]; 
SIG{3} = [SigmaMats.BB(x,1) SigmaMats.BB(x,3); SigmaMats.BB(x,3) SigmaMats.BB(x,2)]; 
MixtureOfGaussiansDraw2dGaussians([MuMats.AA(x,:)' MuMats.AB(x,:)' MuMats.BB(x,:)']' , SIG, labels_vec, legends_vec)
plot(CopyMatA(x), CopyMatB(x), 'm+');
LearnedGenotypesMat(x)

%% [SS II JJ] = intersect(snp_ids(I(Errors)), RLMM.snp_ids(RLMM.SparseSNPs));


figure; hist(max_vec(XXX), 100);

figure; hist(max_vec, 100);


RLMM_Correction = S;


