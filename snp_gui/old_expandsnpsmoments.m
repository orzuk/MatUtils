% 
%
function [mean_mat std_mat std_inv_mat] = ExpandSNPsMoments(RLMM, snp_ids, r_mat);

AssignAllGlobalConstants;


% % Compute also the inverse sigma parameters for the RLMM - this should be done in advance once for the whole RLMM !
% RLMM.SigmaInvMats.AA = RLMM.SigmaMats.AA; RLMM.SigmaInvMats.AB = RLMM.SigmaMats.AB; RLMM.SigmaInvMats.BB = RLMM.SigmaMats.BB;
% for i=1:length(RLMM.snp_ids)
%     M = inv([RLMM.SigmaMats.AA(i,1) RLMM.SigmaMats.AA(i,3); RLMM.SigmaMats.AA(i,3) RLMM.SigmaMats.AA(i,1)]);
%     RLMM.SigmaInvMats.AA(i,1) = M(1,1); RLMM.SigmaInvMats.AA(i,2) = M(2,2); RLMM.SigmaInvMats.AA(i,3) = M(1,2); 
%     M = inv([RLMM.SigmaMats.AB(i,1) RLMM.SigmaMats.AB(i,3); RLMM.SigmaMats.AB(i,3) RLMM.SigmaMats.AA(i,1)]);
%     RLMM.SigmaInvMats.AB(i,1) = M(1,1); RLMM.SigmaInvMats.AB(i,2) = M(2,2); RLMM.SigmaInvMats.AB(i,3) = M(1,2); 
%     M = inv([RLMM.SigmaMats.BB(i,1) RLMM.SigmaMats.BB(i,3); RLMM.SigmaMats.BB(i,3) RLMM.SigmaMats.BB(i,1)]);
%     RLMM.SigmaInvMats.BB(i,1) = M(1,1); RLMM.SigmaInvMats.BB(i,2) = M(2,2); RLMM.SigmaInvMats.BB(i,3) = M(1,2); 
%     if(mod(i,1000) == 0)
%         inverting_snp = i
%     end
% end
 
    
% perform intersection with RLMM to get the correct SNPs moments
[snp_intersect I J] = intersect(snp_ids, RLMM.snp_ids);

num_snps = length(snp_intersect);

% We need to compute the SNP-specific mean and st.d. values
mean_mat = zeros(num_snps, 2*25);
std_mat = ones(num_snps, 3*25);
std_inv_mat = ones(num_snps, 3*25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1: Start with the normal moments
A_copy = 2; B_copy = 0; %AA
mean_mat(I, 1+A_copy+5*B_copy) = RLMM.MuMats.AA(J,1); % A copy
mean_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.MuMats.AA(J,2); % B copy
std_mat(I, 1+A_copy+5*B_copy) = RLMM.SigmaMats.AA(J,1); % A sigma
std_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.SigmaMats.AA(J,2); % B sigma
std_mat(I, 1+A_copy+5*B_copy + 50) = RLMM.SigmaMats.AA(J,3); % covariance
std_inv_mat(I, 1+A_copy+5*B_copy) = RLMM.SigmaInvMats.AA(J,1); % A sigma inv
std_inv_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.SigmaInvMats.AA(J,2); % B sigma inv
std_inv_mat(I, 1+A_copy+5*B_copy + 50) = RLMM.SigmaInvMats.AA(J,3); % covariance inv


A_copy = 1; B_copy = 1; %AB
mean_mat(I, 1+A_copy+5*B_copy) = RLMM.MuMats.AB(J,1); % A copy
mean_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.MuMats.AB(J,2); % B copy
std_mat(I, 1+A_copy+5*B_copy) = RLMM.SigmaMats.AB(J,1); % A sigma
std_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.SigmaMats.AB(J,2); % B sigma
std_mat(I, 1+A_copy+5*B_copy + 50) = RLMM.SigmaMats.AB(J,3); % covariance
std_inv_mat(I, 1+A_copy+5*B_copy) = RLMM.SigmaInvMats.AB(J,1); % A sigma inv
std_inv_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.SigmaInvMats.AB(J,2); % B sigma inv
std_inv_mat(I, 1+A_copy+5*B_copy + 50) = RLMM.SigmaInvMats.AB(J,3); % covariance inv


A_copy = 0; B_copy = 2; %BB
mean_mat(I, 1+A_copy+5*B_copy) = RLMM.MuMats.BB(J,1); % A copy
mean_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.MuMats.BB(J,2); % B copy
std_mat(I, 1+A_copy+5*B_copy) = RLMM.SigmaMats.BB(J,1); % A sigma
std_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.SigmaMats.BB(J,2); % B sigma
std_mat(I, 1+A_copy+5*B_copy + 50) = RLMM.SigmaMats.BB(J,3); % covariance
std_inv_mat(I, 1+A_copy+5*B_copy) = RLMM.SigmaInvMats.BB(J,1); % A sigma inv
std_inv_mat(I, 1+A_copy+5*B_copy + 25) = RLMM.SigmaInvMats.BB(J,2); % B sigma inv
std_inv_mat(I, 1+A_copy+5*B_copy + 50) = RLMM.SigmaInvMats.BB(J,3); % covariance inv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 2: Compute moments for homozygous SNPs
B_copy = 0;
for A_copy = [1 3 4] % A, AAA, AAAA
    mean_mat(:, 1+A_copy+5*B_copy + [0,25]) = r_mat(2,A_copy) * mean_mat(:, 1+2+5*B_copy + [0,25]);
%%    std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = r_mat(2,A_copy) * std_mat(:, 1+2+5*B_copy + [0,25,50]);
%%    std_inv_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = (1/r_mat(2,A_copy)) * std_inv_mat(:, 1+2+5*B_copy + [0,25,50]);
    std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = 1 * std_mat(:, 1+2+5*B_copy + [0,25,50]);
    std_inv_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = 1 * std_inv_mat(:, 1+2+5*B_copy + [0,25,50]);

end
A_copy = 0;
for B_copy = [1 3 4]  % B, BBB, BBBB
    mean_mat(:, 1+A_copy+5*B_copy + [0,25]) = r_mat(2,B_copy) * mean_mat(:, 1+A_copy+5*2 + [0,25]);
%%    std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = r_mat(2,B_copy) * std_mat(:, 1+A_copy+5*2 + [0,25,50]);
%%    std_inv_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = (1/r_mat(2,B_copy)) * std_inv_mat(:, 1+A_copy+5*2 + [0,25,50]);
    std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = 1 * std_mat(:, 1+A_copy+5*2 + [0,25,50]);
    std_inv_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = 1 * std_inv_mat(:, 1+A_copy+5*2 + [0,25,50]);

end

% Special: empty state ! 
A_copy = 0; B_copy = 0; zero_frac = 0.1;  % ratio between the zero copy number and 2 copy number (taken from AB)
mean_mat(:, 1+A_copy+5*B_copy + [0,25]) = zero_frac * mean_mat(:, 1+1+5*1 + [0,25]);
std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = zero_frac * std_mat(:, 1+1+5*1 + [0,25,50]);
std_inv_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = (1/zero_frac) * std_inv_mat(:, 1+1+5*1 + [0,25,50]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 3: Compute moments for heterozygous SNPs
A_copy = 2; B_copy = 2; % AABB
mean_mat(:, 1+A_copy+5*B_copy + [0,25]) = r_mat(2,4) * mean_mat(:, 1+1+5*1 + [0,25]);
%% std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = r_mat(2,4) * std_mat(:, 1+1+5*1 + [0,25,50]);
%% std_inv_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = (1/r_mat(2,4)) * std_inv_mat(:, 1+1+5*1 + [0,25,50]);
std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = 1 * std_mat(:, 1+1+5*1 + [0,25,50]);
std_inv_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = 1 * std_inv_mat(:, 1+1+5*1 + [0,25,50]);


lambda = (r_mat(2,4) - r_mat(2,3)) / (r_mat(2,4)-1); % default lambda (maybe needs to be optimized)  
% lambda = 0; % wrong! remove ! Note that lambda goes to (0,2) and (1-lambda) goes to (1,1)


A_copy = 2; B_copy = 1; % AAB
mean_mat(:, 1+A_copy+5*B_copy + [0,25]) = lambda * mean_mat(:, 1+2+5*0 + [0,25]) + (1-lambda) * r_mat(2,4) * mean_mat(:, 1+1+5*1 + [0,25]);
std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = lambda * std_mat(:, 1+2+5*0 + [0,25,50]) + (1-lambda) * r_mat(2,4) * std_mat(:, 1+1+5*1 + [0,25,50]);
for i=1:num_snps % loop and calculate inverse one by one ..
    tmp_inv = inv( [std_mat(i, 1+A_copy+5*B_copy) std_mat(i, 1+A_copy+5*B_copy+50); ...
        std_mat(i, 1+A_copy+5*B_copy+50) std_mat(i, 1+A_copy+5*B_copy+25)] );
    std_inv_mat(i, 1+A_copy+5*B_copy) = tmp_inv(1,1);
    std_inv_mat(i, 1+A_copy+5*B_copy+25) = tmp_inv(2,2);
    std_inv_mat(i, 1+A_copy+5*B_copy+50) = tmp_inv(1,2);
end
A_copy = 1; B_copy = 2; % ABB
mean_mat(:, 1+A_copy+5*B_copy + [0,25]) = lambda * mean_mat(:, 1+0+5*2 + [0,25]) + (1-lambda) * r_mat(2,4) * mean_mat(:, 1+1+5*1 + [0,25]);
std_mat(:, 1+A_copy+5*B_copy + [0,25,50]) = lambda * std_mat(:, 1+0+5*2 + [0,25,50]) + (1-lambda) * r_mat(2,4) * std_mat(:, 1+1+5*1 + [0,25,50]);
for i=1:num_snps % loop and calculate inverse one by one .. 
    tmp_inv = inv( [std_mat(i, 1+A_copy+5*B_copy) std_mat(i, 1+A_copy+5*B_copy+50); ...
        std_mat(i, 1+A_copy+5*B_copy+50) std_mat(i, 1+A_copy+5*B_copy+25)] );
    std_inv_mat(i, 1+A_copy+5*B_copy) = tmp_inv(1,1);
    std_inv_mat(i, 1+A_copy+5*B_copy+25) = tmp_inv(2,2);
    std_inv_mat(i, 1+A_copy+5*B_copy+50) = tmp_inv(1,2);
end






