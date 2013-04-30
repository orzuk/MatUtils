% Written by Or Zuk 8/2007
%
% Here we compute some parameters which were not present in the original 
% RLMM struct which we learnd from the hapmap samples, such as the inverse
% of the Gaussian covariance matrices, and the determinants.
%
%
% The input:
%
% RLMM - Structure with all parameters
%
% The output:
%
% RLMM - the corrected struct
function RLMM = RLMM_ComputeAuxillaryParams(RLMM)

AssignAllGlobalConstants();

num_snps = length(RLMM.snp_ids);

% Compute also the inverse sigma parameters for the RLMM - this should be done in advance once for the whole RLMM !
RLMM.SigmaInvMats.AA = RLMM.SigmaMats.AA; RLMM.SigmaInvMats.AB = RLMM.SigmaMats.AB; RLMM.SigmaInvMats.BB = RLMM.SigmaMats.BB;
RLMM.SigmaDets.AA = zeros(num_snps, 1); RLMM.SigmaDets.AB = zeros(num_snps, 1); RLMM.SigmaDets.BB = zeros(num_snps, 1);
for i=1:length(RLMM.snp_ids)
    M = inv([RLMM.SigmaMats.AA(i,1) RLMM.SigmaMats.AA(i,3); RLMM.SigmaMats.AA(i,3) RLMM.SigmaMats.AA(i,2)]);
    RLMM.SigmaInvMats.AA(i,1) = M(1,1); RLMM.SigmaInvMats.AA(i,2) = M(2,2); RLMM.SigmaInvMats.AA(i,3) = M(1,2); 
    RLMM.SigmaDets.AA(i) = RLMM.SigmaMats.AA(i,1) * RLMM.SigmaMats.AA(i,2) -  RLMM.SigmaMats.AA(i,3)^2; % get determinant
    M = inv([RLMM.SigmaMats.AB(i,1) RLMM.SigmaMats.AB(i,3); RLMM.SigmaMats.AB(i,3) RLMM.SigmaMats.AB(i,2)]);
    RLMM.SigmaInvMats.AB(i,1) = M(1,1); RLMM.SigmaInvMats.AB(i,2) = M(2,2); RLMM.SigmaInvMats.AB(i,3) = M(1,2); 
    RLMM.SigmaDets.AB(i) = RLMM.SigmaMats.AB(i,1) * RLMM.SigmaMats.AB(i,2) -  RLMM.SigmaMats.AB(i,3)^2; % get determinant
    M = inv([RLMM.SigmaMats.BB(i,1) RLMM.SigmaMats.BB(i,3); RLMM.SigmaMats.BB(i,3) RLMM.SigmaMats.BB(i,2)]);
    RLMM.SigmaInvMats.BB(i,1) = M(1,1); RLMM.SigmaInvMats.BB(i,2) = M(2,2); RLMM.SigmaInvMats.BB(i,3) = M(1,2); 
    RLMM.SigmaDets.BB(i) = RLMM.SigmaMats.BB(i,1) * RLMM.SigmaMats.BB(i,2) -  RLMM.SigmaMats.BB(i,3)^2; % get determinant

    if(mod(i,1000) == 0)
        inverting_snp = i
    end
end


