% Written by Or Zuk 7/2007
%
% Here we implement the RLMM package of Rabbee&Speed (Bioinformatics 2006),
% which was originally written in R, in Matlab.
% This function perform classification - determine the genotypes from the A&B intensitiy
% data. The method is calculating the Mahalahonis distance of each point from 
% all the three Gaussians and taking the shortest one - this is equivalent (up to a small 
% log factor) to taking the maximal-likelihood Gaussian. The function uses the calculated 
% moments of the Multi-Dimensional Gaussian to classify the raw SNP intensities into genotypes.
%
% The input to the function is a test set: matrices of A&B intensities and the calculated moments.
%
% The input:
%
% CopyMatA - Intensities of A allele
% CopyMatB - Intensities of A allele
% MuSNPsMat- The matrix of SNPs mean intensities 
% SigmaSNPsMat- The matrix of SNPs variace intensities 
% MinSparse - The minimum number such that above it the SNP is not considered sparase e.g. 5)
%
% The output:
%
% LearnedGenotypesMat - The inferred genotypes as decided by the classifier
% Mahalahonis - The Mahalahonis distances which were computed and used for classification
function [LearnedGenotypesMat Mahalahonis] = RLMM_ClassifyUsingGaussianParams(CopyMatA, CopyMatB, MuSNPsMat, SigmaSNPsMat);

AssignAllGlobalConstants(); 

NumSNPs = size(CopyMatA, 1);
NumSamples = size(CopyMatA, 2);


% % First loop over all the matrices and compute their inverses
% InvSigmaSNPsMat.AA(1,1,:) = SigmaSNPsMat.AA(:,2);
% InvSigmaSNPsMat.AA(2,2,:) = SigmaSNPsMat.AA(:,1);
% InvSigmaSNPsMat.AA(1,2,:) = -SigmaSNPsMat.AA(:,3);
% InvSigmaSNPsMat.AA(2,1,:) = -SigmaSNPsMat.AA(:,3);
% 
% InvSigmaSNPsMat.AA = InvSigmaSNPsMat.AA ./ ...
%     reshape(repmat(  (SigmaSNPsMat.AA(:,1) .*  SigmaSNPsMat.AA(:,2) - SigmaSNPsMat.AA(:,3).^2), 1, 4)', 2, 2, NumSNPs);  
% InvSigmaSNPsMat.AB(1,1,:) = SigmaSNPsMat.AB(:,2);
% InvSigmaSNPsMat.AB(2,2,:) = SigmaSNPsMat.AB(:,1);
% InvSigmaSNPsMat.AB(1,2,:) = -SigmaSNPsMat.AB(:,3);
% InvSigmaSNPsMat.AB(2,1,:) = -SigmaSNPsMat.AB(:,3);
% InvSigmaSNPsMat.AB = InvSigmaSNPsMat.AB ./ ...
%     reshape(repmat(  (SigmaSNPsMat.AB(:,1) .*  SigmaSNPsMat.AB(:,2) - SigmaSNPsMat.AB(:,3).^2), 1, 4)', 2, 2, NumSNPs);  
% InvSigmaSNPsMat.BB(1,1,:) = SigmaSNPsMat.BB(:,2);
% InvSigmaSNPsMat.BB(2,2,:) = SigmaSNPsMat.BB(:,1);
% InvSigmaSNPsMat.BB(1,2,:) = -SigmaSNPsMat.BB(:,3);
% InvSigmaSNPsMat.BB(2,1,:) = -SigmaSNPsMat.BB(:,3);
% InvSigmaSNPsMat.BB = InvSigmaSNPsMat.BB ./ ...
%     reshape(repmat(  (SigmaSNPsMat.BB(:,1) .*  SigmaSNPsMat.BB(:,2) - SigmaSNPsMat.BB(:,3).^2), 1, 4)', 2, 2, NumSNPs);  


% The above doesn't work well due to singularities. Here we use matlab to
% compute the inverse
InvSigmaSNPsMat.AA(1,1,:) =  SigmaSNPsMat.AA(:,1); InvSigmaSNPsMat.AA(1,2,:) =  SigmaSNPsMat.AA(:,3);
InvSigmaSNPsMat.AA(2,1,:) =  SigmaSNPsMat.AA(:,3); InvSigmaSNPsMat.AA(2,2,:) =  SigmaSNPsMat.AA(:,2);
InvSigmaSNPsMat.AB(1,1,:) =  SigmaSNPsMat.AB(:,1); InvSigmaSNPsMat.AB(1,2,:) =  SigmaSNPsMat.AB(:,3);
InvSigmaSNPsMat.AB(2,1,:) =  SigmaSNPsMat.AB(:,3); InvSigmaSNPsMat.AB(2,2,:) =  SigmaSNPsMat.AB(:,2);
InvSigmaSNPsMat.BB(1,1,:) =  SigmaSNPsMat.BB(:,1); InvSigmaSNPsMat.BB(1,2,:) =  SigmaSNPsMat.BB(:,3);
InvSigmaSNPsMat.BB(2,1,:) =  SigmaSNPsMat.BB(:,3); InvSigmaSNPsMat.BB(2,2,:) =  SigmaSNPsMat.BB(:,2);


for i=1:NumSNPs
    InvSigmaSNPsMat.AA(:,:,i) = inv(InvSigmaSNPsMat.AA(:,:,i)); %   [SigmaSNPsMat.AA(i,1) SigmaSNPsMat.AA(i,3); SigmaSNPsMat.AA(i,3) SigmaSNPsMat.AA(i,2)]);
    InvSigmaSNPsMat.AB(:,:,i) = inv(InvSigmaSNPsMat.AB(:,:,i)); %   [SigmaSNPsMat.AB(i,1) SigmaSNPsMat.AB(i,3); SigmaSNPsMat.AB(i,3) SigmaSNPsMat.AB(i,2)]);
    InvSigmaSNPsMat.BB(:,:,i) = inv(InvSigmaSNPsMat.BB(:,:,i)); %   [SigmaSNPsMat.BB(i,1) SigmaSNPsMat.BB(i,3); SigmaSNPsMat.BB(i,3) SigmaSNPsMat.BB(i,2)]);

    if(mod(i, 1000) == 0)
        do_inv_i_is = i
    end
end

AA_Err = union(find(InvSigmaSNPsMat.AA == -Inf), find(InvSigmaSNPsMat.AA == Inf));
AB_Err = union(find(InvSigmaSNPsMat.AB == -Inf), find(InvSigmaSNPsMat.AB == Inf));
BB_Err = union(find(InvSigmaSNPsMat.BB == -Inf), find(InvSigmaSNPsMat.BB == Inf));

if(~isempty(union(AA_Err, union(AB_Err,BB_Err))))
   Error = 'Inf in Inverse matrices!!!'
   AA_Err
   AB_Err
   BB_Err
end

% Now compute Mahalahonis distances
Mahalahonis.AA = zeros(NumSNPs, NumSamples); 
Mahalahonis.AB = zeros(NumSNPs, NumSamples); 
Mahalahonis.BB = zeros(NumSNPs, NumSamples);

for i=1:NumSNPs
   AA_VEC =  [CopyMatA(i,:)'-MuSNPsMat.AA(i,1) CopyMatB(i,:)'-MuSNPsMat.AA(i,2)];
   AB_VEC =  [CopyMatA(i,:)'-MuSNPsMat.AB(i,1) CopyMatB(i,:)'-MuSNPsMat.AB(i,2)];
   BB_VEC =  [CopyMatA(i,:)'-MuSNPsMat.BB(i,1) CopyMatB(i,:)'-MuSNPsMat.BB(i,2)];

   Mahalahonis.AA(i,:) = sum ( (AA_VEC * InvSigmaSNPsMat.AA(:,:,i)) .* AA_VEC,2 ) + ...
       log(SigmaSNPsMat.AA(i,1)*SigmaSNPsMat.AA(i,2)-SigmaSNPsMat.AA(i,3)*SigmaSNPsMat.AA(i,3)); % add log determinant
   Mahalahonis.AB(i,:) = sum ( (AB_VEC * InvSigmaSNPsMat.AB(:,:,i)) .* AB_VEC,2 ) + ...
    log(SigmaSNPsMat.AB(i,1)*SigmaSNPsMat.AB(i,2)-SigmaSNPsMat.AA(i,3)*SigmaSNPsMat.AB(i,3)); % add log determinant
   Mahalahonis.BB(i,:) = sum ( (BB_VEC * InvSigmaSNPsMat.BB(:,:,i)) .* BB_VEC,2 ) + ...
    log(SigmaSNPsMat.BB(i,1)*SigmaSNPsMat.BB(i,2)-SigmaSNPsMat.BB(i,3)*SigmaSNPsMat.BB(i,3)); % add log determinant
   
   if(mod(i, 1000) == 0)
       i_is = i
   end
end


LearnedGenotypesMat = AA*ones(NumSNPs, NumSamples);
LearnedGenotypesMat(find(Mahalahonis.AB < Mahalahonis.AA)) = AB;
LearnedGenotypesMat(find(Mahalahonis.BB < min(Mahalahonis.AA, Mahalahonis.AB))) = BB;

% add thresholds 
linear_classify_flag=1;
if(linear_classify_flag)
    LearnedGenotypesMat(find(CopyMatB < LINEAR_LOWER_THRESHOLD)) = AA;
    LearnedGenotypesMat(find(CopyMatA < LINEAR_LOWER_THRESHOLD)) = BB;

    LearnedGenotypesMat(find( (CopyMatA > LINEAR_UPPER_THRESHOLD) & (CopyMatB > LINEAR_UPPER_THRESHOLD))) = AB;
end
    
    