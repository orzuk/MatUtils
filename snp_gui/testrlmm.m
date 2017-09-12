% This script is desinged to test the RLMM normalization and SNP genotyping
% functions. 
AssignAllGlobalConstants();

NumSNPs = 500; NumSamples = 1200;

MuMuVec = [2, 0.2, 1, 1, 0.2, 2]; % The six Mu's for AA, AB and BB
% MuMuVec = [12, 0.2, 5, 5, 0.2, 12]; % The six Mu's for AA, AB and BB

pos_flag = 0; pos_ctr=0;
while(~pos_flag)
    pos_ctr=pos_ctr+1
    MuSigmaMat = zeros(6);
%    MuSigmaMat = 0.125 * rand(6); MuSigmaMat = MuSigmaMat + MuSigmaMat';
    MuSigmaMat = 1 .* (MuSigmaMat + diag([0.03, 0.005, 0.015, 0.014, 0.005, 0.03]));

    for i=1:6
        for j=i+1:6
            MuSigmaMat(i,j) = sqrt(MuSigmaMat(i,i) * MuSigmaMat(j,j)) - 0.0000001 *rand(1);
            MuSigmaMat(j,i) = MuSigmaMat(i,j); 
        end
    end
    
    pos_flag = isposdef(MuSigmaMat);
end
% MuSigmaMat = ones(6) + 0.00001 * eye(6);

SigmaMuVec = [0.06 0.03 0.01 0.045 0.045 0.01 0.03 0.06 0.01]; % We need to keep this positive definite
%% SigmaMuVec = [0.09 0.07 0.01 0.08 0.08 0.01 0.07 0.09 0.01]; % We need to keep this positive definite
pos_flag = 0;pos_ctr=0;
while(~pos_flag)
    pos_ctr=pos_ctr+1
    SigmaSigmaMat = 0.0001 * eye(9); % 0.00001 * eye(9);
     for i=1:9
        for j=i+1:9
            SigmaSigmaMat(i,j) = 0.99*sqrt(SigmaSigmaMat(i,i) * SigmaSigmaMat(j,j)) - 0.000000001 *rand(1);
            SigmaSigmaMat(j,i) = SigmaSigmaMat(i,j); 
        end
    end
    
    pos_flag = isposdef(SigmaSigmaMat);
end

[MuSNPsMat, SigmaSNPsMat, CopyMatA, CopyMatB, GenotypesMat] = ...
    SimulateSNPsMats(NumSNPs, NumSamples, MuMuVec, MuSigmaMat,SigmaMuVec, SigmaSigmaMat);

% Call the normalization functions: 
PsuedoCount=40;
RLMM = RLMM_LearnGaussianParams(CopyMatA, CopyMatB, GenotypesMat, NumSamples*0.95/3, PsuedoCount);
LearnedMuSNPsMat = RLMM.MuMats 
LearnedSigmaSNPsMat = RLMM.SigmaMats 
SparseSNPs = RLMM.SparseSNPs;


% Call the genotype Classification functions with the learned parameters
% (should be on a test set which is different from the training set used
% for the parameter estimation) 

%[InferredGenotypesMat] = RLMM_ClassifyUsingGaussianParams(CopyMatA, CopyMatB, MuSNPsMat, SigmaSNPsMat);

[LearnedGenotypesMat Mahalahonis] = RLMM_ClassifyUsingGaussianParams(CopyMatA, CopyMatB, LearnedMuSNPsMat, LearnedSigmaSNPsMat);

% Compare the learned Mew_mats to the original mats
figure; plot(MuSNPsMat.AA(:), LearnedMuSNPsMat.AA(:), '.'); 
title('True \mu vs. Learned \mu values'); xlabel('True'); ylabel('Learned');

figure; hold on; plot(reshape(SigmaSNPsMat.AA(1,1,:), 1, NumSNPs), LearnedSigmaSNPsMat.AA(:,1), '.'); 
plot(reshape(SigmaSNPsMat.AA(2,2,:), 1, NumSNPs), LearnedSigmaSNPsMat.AA(:,2), 'r.'); 
plot(reshape(SigmaSNPsMat.AA(1,2,:), 1, NumSNPs), LearnedSigmaSNPsMat.AA(:,3), 'g.');
title('AA True \Sigma vs. Learned \Sigma values'); xlabel('True'); ylabel('Learned');


figure; hold on; plot(reshape(SigmaSNPsMat.AB(1,1,:), 1, NumSNPs), LearnedSigmaSNPsMat.AB(:,1), '.'); 
plot(reshape(SigmaSNPsMat.AB(2,2,:), 1, NumSNPs), LearnedSigmaSNPsMat.AB(:,2), 'r.'); 
plot(reshape(SigmaSNPsMat.AB(1,2,:), 1, NumSNPs), LearnedSigmaSNPsMat.AB(:,3), 'g.');
title('AB True \Sigma vs. Learned \Sigma values'); xlabel('True'); ylabel('Learned');

figure; hold on; plot(reshape(SigmaSNPsMat.BB(1,1,:), 1, NumSNPs), LearnedSigmaSNPsMat.BB(:,1), '.'); 
plot(reshape(SigmaSNPsMat.BB(2,2,:), 1, NumSNPs), LearnedSigmaSNPsMat.BB(:,2), 'r.'); 
plot(reshape(SigmaSNPsMat.BB(1,2,:), 1, NumSNPs), LearnedSigmaSNPsMat.BB(:,3), 'g.');
title('BB True \Sigma vs. Learned \Sigma values'); xlabel('True'); ylabel('Learned');
figure; plot(reshape(SigmaSNPsMat.AA(1,2,:), 1, NumSNPs), LearnedSigmaSNPsMat.AA(:,3), '.'); 
title('True COV vs. Learned COV values'); xlabel('True'); ylabel('Learned');

% look at the data from one SNP
snp_ind = 62;
labels_vec = {'A Intensity', 'B Intensity'}; legends_vec = {'AA', 'AB', 'BB'};
FirstMuVecs = [LearnedMuSNPsMat.AA(snp_ind,:)' LearnedMuSNPsMat.AB(snp_ind,:)' LearnedMuSNPsMat.BB(snp_ind,:)']'; 
FirstSigmaMats{1} = [LearnedSigmaSNPsMat.AA(snp_ind,1) LearnedSigmaSNPsMat.AA(snp_ind,3); ...
    LearnedSigmaSNPsMat.AA(snp_ind,3) LearnedSigmaSNPsMat.AA(snp_ind,2)];
FirstSigmaMats{2} = [LearnedSigmaSNPsMat.AB(snp_ind,1) LearnedSigmaSNPsMat.AB(snp_ind,3); ...
    LearnedSigmaSNPsMat.AB(snp_ind,3) LearnedSigmaSNPsMat.AB(snp_ind,2)];
FirstSigmaMats{3} = [LearnedSigmaSNPsMat.BB(snp_ind,1) LearnedSigmaSNPsMat.BB(snp_ind,3); ...
    LearnedSigmaSNPsMat.BB(snp_ind,3) LearnedSigmaSNPsMat.BB(snp_ind,2)];
MixtureOfGaussiansDraw2dGaussians(FirstMuVecs, FirstSigmaMats, labels_vec, legends_vec);
PlotAlleleRatios(CopyMatA(snp_ind,:)+CopyMatB(snp_ind,:), CopyMatA(snp_ind,:)./CopyMatB(snp_ind,:), ...
    GenotypesMat(snp_ind,:), ['True  SNP ' num2str(snp_ind) ' Genotypes'])
PlotAlleleRatios(CopyMatA(snp_ind,:)+CopyMatB(snp_ind,:), CopyMatA(snp_ind,:)./CopyMatB(snp_ind,:), ...
    LearnedGenotypesMat(snp_ind,:), ['Learned  SNP ' num2str(snp_ind) ' Genotypes'])

NumSparseSNPs = length(SparseSNPs)
ClassificationError = sum(sum(LearnedGenotypesMat ~= GenotypesMat)) ./ (NumSNPs * NumSamples)

DenseSNPs = setdiff(1:NumSNPs,SparseSNPs);
ClassificationErrorDense = sum(sum(LearnedGenotypesMat(DenseSNPs,:) ~= GenotypesMat(DenseSNPs,:))) ./ (length(DenseSNPs) * NumSamples)
ClassificationErrorSparse = sum(sum(LearnedGenotypesMat(SparseSNPs,:) ~= GenotypesMat(SparseSNPs,:))) ./ (length(SparseSNPs) * NumSamples)

