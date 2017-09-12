% Auxillary function. Simuale data from multi-variate Gaussian
% distribution, where for each SNP we've got different Gaussians parameters
% for generating it's intensity data. We also randomize the genotypes
% (currently completely at random without any structure) 
%
function [MuSNPsMat, SigmaSNPsMat, CopyMatA, CopyMatB, GenotypesMat] = SimulateSNPsMats(NumSNPs, NumSamples, ...
    MuMuVec, MuSigmaMat,SigmaMuVec, SigmaSigmaMat)

AssignAllGlobalConstants();  EPSILON = 0.0000000001;

% Randomize SNP Intensity moments
isposdef(MuSigmaMat)
isposdef(SigmaSigmaMat)
TmpMuSNPsMat = mgd(NumSNPs,6,MuMuVec,MuSigmaMat); 
TmpSigmaSNPsMat = max(EPSILON, mgd(NumSNPs,9,SigmaMuVec,SigmaSigmaMat));

MuSNPsMat.AA = TmpMuSNPsMat(:,1:2);
MuSNPsMat.AB = TmpMuSNPsMat(:,3:4);
MuSNPsMat.BB = TmpMuSNPsMat(:,5:6);

SigmaSNPsMat.AA(1,:,:) = TmpSigmaSNPsMat(:,[1,3])';
SigmaSNPsMat.AA(2,:,:) = TmpSigmaSNPsMat(:,[3,2])';
SigmaSNPsMat.AB(1,:,:) = TmpSigmaSNPsMat(:,[4,6])';
SigmaSNPsMat.AB(2,:,:) = TmpSigmaSNPsMat(:,[6,5])';
SigmaSNPsMat.BB(1,:,:) = TmpSigmaSNPsMat(:,[4,6])';
SigmaSNPsMat.BB(2,:,:) = TmpSigmaSNPsMat(:,[6,5])';


% Randomize Genotypes
GenotypesMat = ceil(3*rand(NumSNPs, NumSamples));
AA_Counts =  sum(GenotypesMat == AA,2);
AB_Counts =  sum(GenotypesMat == AB,2);
BB_Counts =  sum(GenotypesMat == BB,2);

AA_CountMat = zeros(NumSNPs, NumSamples); AA_CountMat(find(GenotypesMat == AA)) = 1;
AB_CountMat = zeros(NumSNPs, NumSamples); AB_CountMat(find(GenotypesMat == AB)) = 1;
BB_CountMat = zeros(NumSNPs, NumSamples); BB_CountMat(find(GenotypesMat == BB)) = 1;


% Randomize intensities - do it in a primitive looping way
CopyMatA = zeros(NumSamples, NumSNPs); CopyMatB = zeros(NumSamples, NumSNPs);

for i=1:NumSNPs
    if(~isposdef(SigmaSNPsMat.AA(:,:,i)))
        SigmaSNPsMat.AA(1,1,i) = abs( SigmaSNPsMat.AA(1,1,i));
        SigmaSNPsMat.AA(2,2,i) =  2 * SigmaSNPsMat.AA(1,2,i)^2 / SigmaSNPsMat.AA(1,1,i);
        NOTPOS = i %       break;
    end
    if(~isposdef(SigmaSNPsMat.AB(:,:,i)))
        SigmaSNPsMat.AB(1,1,i) = abs( SigmaSNPsMat.AB(1,1,i));
        SigmaSNPsMat.AB(2,2,i) =  2 * SigmaSNPsMat.AB(1,2,i)^2 / SigmaSNPsMat.AB(1,1,i);
        NOTPOS = i %       break;
    end
    if(~isposdef(SigmaSNPsMat.BB(:,:,i)))
        SigmaSNPsMat.BB(1,1,i) = abs( SigmaSNPsMat.BB(1,1,i));
        SigmaSNPsMat.BB(2,2,i) =  2 * SigmaSNPsMat.BB(1,2,i)^2 / SigmaSNPsMat.BB(1,1,i);
        NOTPOS = i    %    break;
    end    
    TmpMatAA = mgd(AA_Counts(i,:),2, MuSNPsMat.AA(i,:), SigmaSNPsMat.AA(:,:,i));
    TmpMatAB = mgd(AB_Counts(i,:),2, MuSNPsMat.AB(i,:), SigmaSNPsMat.AB(:,:,i));
    TmpMatBB = mgd(BB_Counts(i,:),2, MuSNPsMat.BB(i,:), SigmaSNPsMat.BB(:,:,i));

    if( (~isreal(TmpMatAA)) | (~isreal(TmpMatAB)) | (~isreal(TmpMatBB)) )
        break;
    end
    
    CopyMatA(find(AA_CountMat(i,:)),i) = TmpMatAA(:,1); 
    CopyMatB(find(AA_CountMat(i,:)),i) = TmpMatAA(:,2); 
    CopyMatA(find(AB_CountMat(i,:)),i) = TmpMatAB(:,1); 
    CopyMatB(find(AB_CountMat(i,:)),i) = TmpMatAB(:,2); 
    CopyMatA(find(BB_CountMat(i,:)),i) = TmpMatBB(:,1); 
    CopyMatB(find(BB_CountMat(i,:)),i) = TmpMatBB(:,2); 

    if(mod(i,100) == 0)
        i_is = i
    end
    
end

CopyMatA = CopyMatA'; CopyMatB = CopyMatB';
