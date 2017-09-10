% Written by Or Zuk 5/2007
%
% This function computes the most probable values of alpha&beta
% copy-number, genotypes etc. from the Gamma probs. this is returned in VitStruct.
% The function also computes the marginal probabilities of alpha&beta
% copynumber, genotype etc. This is returned in ProbsStruct.
% The method is simply to sum-up all the relevant probs and take the
% maximal one whenever neccessary.
% The convention for the probability (Gamma_Probs) table is as follows :
% __________________________________________________________
% | alpha-copy | alpha-genotype | beta-copy | beta-genotype |
% |_________________________________________________________|
% So the lsb is beta-genotypes (1 bit), then 3 possibilities for beta-copy,
% then alpha-genotype (1 bit) and finally alpha-copy (3 possibilities)
%
% Input:
% Gamma_Probs - A table of 36*seq_len of states probs
% HMM_MODEL - Structure containing all model's parameters
% do_couples - a flag which says if we want to consider also correlations between
% adjacent pairs (its default value should be one - much slower but more accurate and less 'bumpy')
% joint_flag - If it's on, we take the A&B and alpha&beta values from the joint most probable
% copy # values, which is recommended since otherwise we may get incorrect
% copy # by small random fluctuations. Note: The function works poorly when this flag is off!!!
%
%
% Output :
% VitStruct - A structure containing the clumped (discrete) values for copy
% # and genotypes (seperately for each chromosome and jointly)
% ProbsStruct - Table of probabilities projected on genotype/copy #
%
function [VitStruct ProbsStruct] = GetBestMarginalPredictions(Gamma_Probs, HMM_MODEL, do_couples, joint_flag)

x_dim = HMM_MODEL.x_dim; x_dim_eq = 2*x_dim-1;

% Tables for indexing the A and B copy
A_copy_tab = [ 0, 0, 1, 0, 2, 0, 0, 0, 1, ...
    0, 2, 0, 1, 1, 2, 1, 3, 1, ...
    0, 0, 1, 0, 2, 0, 2, 2, 3, ...
    2, 4, 2, 0, 0, 1, 0, 2, 0];
B_copy_tab = [ 0, 0, 0, 1, 0, 2, 0, 0, 0, ...
    1, 0, 2, 0, 0, 0, 1, 0, 2, ...
    1, 1, 1, 2, 1, 3, 0, 0, 0, ...
    1, 0, 2, 2, 2, 2, 3, 2, 4];
AB_copy_tab = A_copy_tab + x_dim_eq*B_copy_tab; % here put A in the lsb (base 5) and B in the  msb (another 5)

% Check if Gamma is 'compressed' and un-compress it. First calculate the probs. % New!! No need for compression !!!
RelevantInds = [1:6 17:22 33:38 49:54 65:70 81:86];
if(size(Gamma_Probs,1) == 36)
    TempGammaProbs = Gamma_Probs;
else
    TempGammaProbs = Gamma_Probs(RelevantInds,:);
end

TOL = 0.00000000000000000000000000000000000000000001;
% First get the participating indexes:
n = size(Gamma_Probs,2);


ProbsStruct.alpha_genotype = [sum(TempGammaProbs(1:2:35,:))' sum(TempGammaProbs(2:2:36,:))']; % alpha genotypes are lsb
ProbsStruct.beta_genotype = [sum(TempGammaProbs([1:6 13:18 25:30],:))' sum(TempGammaProbs([7:12 19:24 31:36],:))']; % beta genotypes are 3rd
ProbsStruct.joint_genotype = [sum(TempGammaProbs([1:2:5 13:2:17 25:2:29],:))' ...
    sum(TempGammaProbs([2:2:6 14:2:18 26:2:30],:))' ...
    sum(TempGammaProbs([7:2:11 19:2:23 31:2:35],:))' ...
    sum(TempGammaProbs([8:2:12 20:2:24 32:2:36],:))'];
ProbsStruct.alpha_copy = [sum(TempGammaProbs([1:2 7:8 13:14 19:20 25:26 31:32],:))' ...
    sum(TempGammaProbs([3:4 9:10 15:16 21:22 27:28 33:34],:))' ...
    sum(TempGammaProbs([5:6 11:12 17:18 23:24 29:30 35:36],:))']; % alpha copy are 2nd
ProbsStruct.beta_copy = [sum(TempGammaProbs([1:12],:))' sum(TempGammaProbs([13:24],:))' sum(TempGammaProbs([25:36],:))']; % beta copy are msb
ProbsStruct.joint_copy = TempGammaProbs([1:2:5 13:2:17 25:2:29],:) + ...
    TempGammaProbs([2:2:6 14:2:18 26:2:30],:) + ...
    TempGammaProbs([7:2:11 19:2:23 31:2:35],:) + ...
    TempGammaProbs([8:2:12 20:2:24 32:2:36],:);
ProbsStruct.total_copy = [ProbsStruct.joint_copy(1,:)' sum(ProbsStruct.joint_copy([2,4],:))' ...
    sum(ProbsStruct.joint_copy([3,5,7],:))' sum(ProbsStruct.joint_copy([6,8],:))' ProbsStruct.joint_copy(9,:)'];
ProbsStruct.joint_copy = ProbsStruct.joint_copy';
% New: Add A&B copy - how ?
ProbsStruct.A_copy = zeros(n,x_dim_eq); ProbsStruct.B_copy = zeros(n,x_dim_eq); ProbsStruct.AB_copy = zeros(n,x_dim_eq^2);
for i=1:x_dim_eq
    A_Inds = find(A_copy_tab == i-1);  B_Inds = find(B_copy_tab == i-1);
    ProbsStruct.A_copy(:,i) = sum(TempGammaProbs(A_Inds,:), 1);
    ProbsStruct.B_copy(:,i) = sum(TempGammaProbs(B_Inds,:), 1);
    for j=1:x_dim_eq
        B_Inds = find(B_copy_tab == j-1); AB_Inds = intersect(A_Inds, B_Inds);
        if(~isempty(AB_Inds))
            ProbsStruct.AB_copy(:,i+x_dim_eq*(j-1)+1) =  sum(TempGammaProbs(AB_Inds,:), 1);
        end
    end
end

if(do_couples == 0) % Perform the easy marginal inference
    V_Inds = zeros(1,2*2*x_dim*x_dim);
    VitStruct.alpha_copy = zeros(x_dim, n);
    VitStruct.beta_copy = zeros(x_dim, n);
    VitStruct.joint_copy = zeros(x_dim*x_dim, n); % Why not just zeros(1,n)
    VitStruct.alpha_genotype = zeros(2, n);
    VitStruct.beta_genotype = zeros(2, n);
    VitStruct.joint_genotype = zeros(2*2, n);
    VitStruct.A_copy = zeros(x_dim_eq, n);
    VitStruct.B_copy = zeros(x_dim_eq, n);
    VitStruct.AB_copy = zeros( x_dim_eq^2, n);

    i=1;

    % Loop in the correct order
    for j1 = 0:x_dim-1
        for j_geno1 = 0:1
            for j2 = 0:x_dim-1
                for j_geno2 = 0:1
                    V_Inds(i) = j_geno1 + 2*j1 + 6*j_geno2 + 12*j2 + 1; % The extra 1 is since we srart with zero and not one

                    VitStruct.alpha_copy(j1+1,:) =  VitStruct.alpha_copy(j1+1,:) + TempGammaProbs(V_Inds(i),:); % update alpha copy number
                    VitStruct.beta_copy(j2+1,:) =  VitStruct.beta_copy(j2+1,:) + TempGammaProbs(V_Inds(i),:); % update beta copy number
                    VitStruct.joint_copy(j1+j2*x_dim+1,:) =  VitStruct.joint_copy(j1+j2*x_dim+1,:) + TempGammaProbs(V_Inds(i),:); % update alpha copy number
                    VitStruct.alpha_genotype(j_geno1+1,:) =  VitStruct.alpha_genotype(j_geno1+1,:) + TempGammaProbs(V_Inds(i),:); % update alpha genotype
                    VitStruct.beta_genotype(j_geno2+1,:) =  VitStruct.beta_genotype(j_geno2+1,:) + TempGammaProbs(V_Inds(i),:); % update beta genotype
                    VitStruct.joint_genotype(j_geno1+2*j_geno2+1,:) = VitStruct.joint_genotype(j_geno1+2*j_geno2+1,:) + TempGammaProbs(V_Inds(i),:); % update beta genotype

                    VitStruct.A_copy(A_copy_tab(V_Inds(i))+1,:) =  VitStruct.A_copy(A_copy_tab(V_Inds(i))+1,:) + TempGammaProbs(V_Inds(i),:); % update alpha copy number
                    VitStruct.B_copy(B_copy_tab(V_Inds(i))+1,:) =  VitStruct.B_copy(B_copy_tab(V_Inds(i))+1,:) + TempGammaProbs(V_Inds(i),:); % update alpha copy number
                    VitStruct.AB_copy(AB_copy_tab(V_Inds(i))+1,:) =  VitStruct.AB_copy(AB_copy_tab(V_Inds(i))+1,:) + TempGammaProbs(V_Inds(i),:); % update alpha copy number

                    i=i+1;
                end
            end
        end
    end
    % Now get the alpha and beta genotypes, simply by taking the maximum
    [Dummy VitStruct.alpha_copy] = max(VitStruct.alpha_copy); VitStruct.alpha_copy = (VitStruct.alpha_copy - 1)';
    [Dummy VitStruct.beta_copy] = max(VitStruct.beta_copy); VitStruct.beta_copy = (VitStruct.beta_copy - 1)';
    [Dummy VitStruct.alpha_genotype] = max(VitStruct.alpha_genotype); VitStruct.alpha_genotype = (VitStruct.alpha_genotype - 1)';
    [Dummy VitStruct.beta_genotype] = max(VitStruct.beta_genotype); VitStruct.beta_genotype = (VitStruct.beta_genotype - 1)';
    [Dummy VitStruct.joint_genotype] = max(VitStruct.joint_genotype); VitStruct.joint_genotype = (VitStruct.joint_genotype - 1)';
    [Dummy VitStruct.joint_copy] = max(VitStruct.joint_copy); VitStruct.joint_copy = (VitStruct.joint_copy - 1)';     % Now try to get alpha&beta together. Should we leave this as is? or take alpha and beta seperately?

    [Dummy VitStruct.A_copy] = max(VitStruct.A_copy); VitStruct.A_copy = (VitStruct.A_copy - 1)';
    [Dummy VitStruct.B_copy] = max(VitStruct.B_copy); VitStruct.B_copy = (VitStruct.B_copy - 1)';
    [Dummy VitStruct.AB_copy] = max(VitStruct.AB_copy); VitStruct.AB_copy = (VitStruct.AB_copy - 1)';

    % Here it is assumed that genotype 0 is A and 1 is B - this is the
    % corrected convention! We take everything from the joint over the two chromosomes
    if(joint_flag)
        % Now try to get alpha&beta together. Should we leave this as is? or take alpha and beta seperately?
        VitStruct.alpha_copy = mod(VitStruct.joint_copy, x_dim);
        VitStruct.beta_copy = (VitStruct.joint_copy-VitStruct.alpha_copy)/x_dim;
        VitStruct.alpha_genotype = bitget(VitStruct.joint_genotype,1);
        VitStruct.beta_genotype = bitget(VitStruct.joint_genotype,2);
        VitStruct.A_copy = mod(VitStruct.AB_copy, x_dim_eq); % A is the lsb
        VitStruct.B_copy = floor(VitStruct.AB_copy/x_dim_eq); % B is the msb
    end

    VitStruct.total_copy = VitStruct.alpha_copy + VitStruct.beta_copy;
    VitStruct.total_copy_from_AB = VitStruct.A_copy + VitStruct.B_copy;

else  % Do a more sophisticated inference based on adjacent couples
    V_Inds = zeros(1,2*2*x_dim*x_dim);
    V_Inds2 = zeros(1,2*2*x_dim*x_dim);
    VitStruct.alpha_copy = zeros(x_dim,x_dim, n);
    VitStruct.beta_copy = zeros(x_dim,x_dim, n);
%%    VitStruct.joint_copy_couples = zeros(x_dim^2, x_dim^2, n);
%%    VitStruct.AB_copy_couples = zeros(x_dim_eq^2, x_dim_eq^2, n);
    VitStruct.alpha_genotype = zeros(2, n);
    VitStruct.beta_genotype = zeros(2, n);
    VitStruct.joint_genotype = zeros(2*2, n);

    i=1;
    for j_geno1 = 0:1
        for j_geno2 = 0:1
            for j1 = 0:x_dim-1
                for j2 = 0:x_dim-1
                    j=1;
                    V_Inds(i) = j_geno1 + 2*j1 + 6*j_geno2 + 12*j2 + 1; % The extra 1 is since we srart with zero and not one
                    VitStruct.joint_genotype(j_geno1+2*j_geno2+1,:) = VitStruct.joint_genotype(j_geno1+2*j_geno2+1,:) + ...
                        TempGammaProbs(V_Inds(i),:); % update joint genotype
                    % % %                     for k_geno1 = 0:1
                    % % %                         for k_geno2 = 0:1
                    % % %                             for k1 = 0:x_dim-1
                    % % %                                 for k2 = 0:x_dim-1
                    % % %                                     V_Inds(i) = j_geno1 + 2*j1 + 6*j_geno2 + 12*j2 + 1; % The extra 1 is since we srart with zero and not one
                    % % %                                     V_Inds2(j) = k_geno1 + 2*k1 + 6*k_geno2 + 12*k2 + 1; % The extra 1 is since we srart with zero and not one
                    % % %                                     VitStruct.joint_copy_couples(j1+j2*x_dim+1,k1+k2*x_dim+1,1:end-1) = ...
                    % % %                                         VitStruct.joint_copy_couples(j1+j2*x_dim+1,k1+k2*x_dim+1,1:end-1) + ...
                    % % %                                         reshape(TempGammaProbs(V_Inds(i),1:end-1) .* TempGammaProbs(V_Inds2(j),2:end) .* ...
                    % % %                                         HMM_MODEL.M(j1+1,k1+1) .* HMM_MODEL.M(j2+1,k2+1) .* ...
                    % % %                                         HMM_MODEL.PLACE_M(1:n-1,j_geno1+2*k_geno1+1)' .* ...
                    % % %                                         HMM_MODEL.PLACE_M(1:n-1,j_geno2+2*k_geno2+1)', 1, 1, n-1); % update alpha-beta copy number. Last two rows represent the interaction term
                    % % %                                     VitStruct.joint_copy_couples(j1+j2*x_dim+1,k1+k2*x_dim+1,end) = ...
                    % % %                                         VitStruct.joint_copy_couples(j1+j2*x_dim+1,k1+k2*x_dim+1,end) + ...
                    % % %                                         TempGammaProbs(V_Inds2(j),end);  % Update the last according to the marginal
                    % % %
                    % % %                                     VitStruct.AB_copy_couples(AB_copy_tab(V_Inds(i))+1,AB_copy_tab(V_Inds2(j))+1,1:end-1) = ...
                    % % %                                         VitStruct.AB_copy_couples(AB_copy_tab(V_Inds(i))+1,AB_copy_tab(V_Inds2(j))+1,1:end-1) + ...
                    % % %                                         reshape(TempGammaProbs(V_Inds(i),1:end-1) .* TempGammaProbs(V_Inds2(j),2:end) .* ...
                    % % %                                         HMM_MODEL.M(j1+1,k1+1) .* HMM_MODEL.M(j2+1,k2+1) .* ...
                    % % %                                         HMM_MODEL.PLACE_M(1:n-1,j_geno1+2*k_geno1+1)' .* ...
                    % % %                                         HMM_MODEL.PLACE_M(1:n-1,j_geno2+2*k_geno2+1)', 1, 1, n-1); % update AB copy number. Last two rows represent the interaction term
                    % % %                                     VitStruct.AB_copy_couples(AB_copy_tab(V_Inds(i))+1,AB_copy_tab(V_Inds2(j))+1,end) = ...
                    % % %                                         VitStruct.AB_copy_couples(AB_copy_tab(V_Inds(i))+1,AB_copy_tab(V_Inds2(j))+1,end) + ...
                    % % %                                         TempGammaProbs(V_Inds2(j),end);  % Update the last according to the marginal
                    % % %
                    % % %                                     j=j+1;
                    % % %                                 end
                    % % %                             end
                    % % %                         end
                    % % %                     end
                    i=i+1;
                end
            end
        end
    end
    % % %     [Dummy VitStruct.alpha_genotype] = max(VitStruct.alpha_genotype); VitStruct.alpha_genotype = (VitStruct.alpha_genotype - 1)';
    % % %     [Dummy VitStruct.beta_genotype] = max(VitStruct.beta_genotype); VitStruct.beta_genotype = (VitStruct.beta_genotype - 1)';
    [Dummy VitStruct.joint_genotype] = max(VitStruct.joint_genotype); VitStruct.joint_genotype = (VitStruct.joint_genotype - 1)';
    % % %
    % % %     % Now try to get alpha&beta together
    % % %     % Normalize: Transfer joint probs into conditional probs:
    % % %     % % %     VitStruct.joint_copy_couples(:,:,end) = 1/(x_dim^4);
    % % %     VitStruct.joint_copy = sum(VitStruct.joint_copy_couples, 2);
    % % %     for j=1:x_dim^2
    % % %         VitStruct.joint_copy_couples(:,j,:) = (VitStruct.joint_copy_couples(:,j,:)+TOL) ./ (VitStruct.joint_copy+x_dim^2*TOL);
    % % %     end
    % % %     VitStruct.AB_copy = sum(VitStruct.AB_copy_couples, 2);
    % % %     for j=1:x_dim_eq^2
    % % %         VitStruct.AB_copy_couples(:,j,:) = (VitStruct.AB_copy_couples(:,j,:)+TOL) ./ (VitStruct.AB_copy+x_dim_eq^2*TOL);
    % % %     end
    % % %
    % % %
    % % %     % Do a 'Viterbi' algorithm again!
    % % %     delta = zeros(x_dim^2,n); phi = zeros(x_dim^2, n); VitStruct.joint_copy2 = zeros(1,n);
    % % %     delta_AB = zeros(x_dim_eq^2,n); phi_AB = zeros(x_dim_eq^2,n);  VitStruct.AB_copy2 = zeros(1,n);
    % % %     % Initilization
    % % %     delta(:,1) = VitStruct.joint_copy(:,1);  phi(:,1) = -1; % Should be irrelevant!
    % % %     delta_AB(:,1) = VitStruct.AB_copy(:,1);  phi_AB(:,1) = -1; % Should be irrelevant!
    % % %     % Recursion
    % % %     for t=2:n
    % % %         for j=1:x_dim^2
    % % %             [delta(j, t) phi(j,t) ] = max(delta(:,t-1) .* ...
    % % %                 reshape(VitStruct.joint_copy_couples(:,j,t-1), x_dim^2, 1)); % Keep the delta value and the path
    % % %         end
    % % %         delta(:,t) = delta(:,t) ./ sum(delta(:,t));         % Do scaling
    % % %
    % % %         for j=1:x_dim_eq^2
    % % %             [delta_AB(j, t) phi_AB(j,t) ] = max(delta_AB(:,t-1) .* ...
    % % %                 reshape(VitStruct.AB_copy_couples(:,j,t-1), x_dim_eq^2, 1)); % Keep the delta_AB value and the path
    % % %         end
    % % %         delta_AB(:,t) = delta_AB(:,t) ./ sum(delta_AB(:,t));         % Do scaling for AB
    % % %     end
    % % %     phi = phi-1; phi_AB = phi_AB-1;
    % % %     % Termination
    % % %     [delta_star q_star] = max(delta(:,n)); VitStruct.joint_copy2(end) = q_star-1; % try ..
    % % %     [delta_star_AB q_star_AB] = max(delta_AB(:,n)); VitStruct.AB_copy2(end) = q_star_AB-1; % try AB ..
    % % %     % Backtracking
    % % %     for t=n-1:-1:1
    % % %         VitStruct.joint_copy2(t) = phi(VitStruct.joint_copy2(t+1)+1,t+1);
    % % %         VitStruct.AB_copy2(t) = phi_AB(VitStruct.AB_copy2(t+1)+1,t+1);
    % % %     end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New: Simpler Approach: %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare Inds tables
    IndsTab{1} = [1 2 7 8]; IndsTab{2} = [3 4 9 10]; IndsTab{3} = [5 6 11 12];
    IndsTab{4} = [13 14 19 20]; IndsTab{5} = [15 16 21 22]; IndsTab{6} = [17 18 23 24];
    IndsTab{7} = [25 26 31 32]; IndsTab{8} = [27 28 33 34]; IndsTab{9} = [29 30 35 36];
    for i=1:36
        IndsTabAB{i} = find(AB_copy_tab == i-1);
    end
    NonEmptyInds = unique(AB_copy_tab)+1;
    ShortGammaVec = zeros(1,36);

    TempGammaProbs2 = TempGammaProbs;
    delta2 = zeros(1,n);
    [Dummy delta2(1)] = max(TempGammaProbs([1:2:5 13:2:17 25:2:29],1) + ...
        TempGammaProbs([2:2:6 14:2:18 26:2:30],1) + ...
        TempGammaProbs([7:2:11 19:2:23 31:2:35],1) + ...
        TempGammaProbs([8:2:12 20:2:24 32:2:36],1));
    for t=2:n
        GammaProbsSum = sum( TempGammaProbs(IndsTab{delta2(t-1)},t) );
        TempGammaProbs2(:,t) = TempGammaProbs2(:,t) ./ GammaProbsSum;
        for j=1:36
            j_geno1 = mod(j,2); j_geno2 = mod(floor(j/6),2); j1 = mod(floor(j/2),3);  j2 = mod(floor(j/12),3);
            ProdTerm=0;
            for i=IndsTab{delta2(t-1)}
                i_geno1 = mod(i,2); i_geno2 = mod(floor(i/6),2); i1 = mod(floor(i/2),3);  i2 = mod(floor(i/12),3);
                ProdTerm = ProdTerm + TempGammaProbs(i,t) * ...
                    HMM_MODEL.M(i1+1,j1+1) .* HMM_MODEL.M(i2+1,j2+1) .* ...
                    HMM_MODEL.PLACE_M(t-1,i_geno1+2*j_geno1+1) .* ...
                    HMM_MODEL.PLACE_M(t-1,i_geno2+2*j_geno2+1);
            end
            TempGammaProbs2(j,t) = TempGammaProbs2(j,t) .* ProdTerm;
        end
        [Dummy delta2(t)] = max(TempGammaProbs2([1:2:5 13:2:17 25:2:29],t) + ...
            TempGammaProbs2([2:2:6 14:2:18 26:2:30],t) + ...
            TempGammaProbs2([7:2:11 19:2:23 31:2:35],t) + ...
            TempGammaProbs2([8:2:12 20:2:24 32:2:36],t));
    end
    VitStruct.joint_copy2 = delta2-1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TempGammaProbs2 = TempGammaProbs;
    delta2 = zeros(1,n);
    for i=NonEmptyInds
        ShortGammaVec(i) = sum(TempGammaProbs(IndsTabAB{i},1));
    end
    [Dummy delta2(1)] = max(ShortGammaVec);
    for t=2:n
        GammaProbsSum = sum( TempGammaProbs(IndsTabAB{delta2(t-1)},t) );
        TempGammaProbs2(:,t) = TempGammaProbs2(:,t) ./ GammaProbsSum;
        for j=1:36
            j_geno1 = mod(j,2); j_geno2 = mod(floor(j/6),2); j1 = mod(floor(j/2),3);  j2 = mod(floor(j/12),3);
            ProdTerm=0;
            for i=IndsTabAB{delta2(t-1)}
                i_geno1 = mod(i,2); i_geno2 = mod(floor(i/6),2); i1 = mod(floor(i/2),3);  i2 = mod(floor(i/12),3);
                ProdTerm = ProdTerm + TempGammaProbs(i,t) * ...
                    HMM_MODEL.M(i1+1,j1+1) .* HMM_MODEL.M(i2+1,j2+1) .* ...
                    HMM_MODEL.PLACE_M(t-1,i_geno1+2*j_geno1+1) .* ...
                    HMM_MODEL.PLACE_M(t-1,i_geno2+2*j_geno2+1);
            end
            TempGammaProbs2(j,t) = TempGammaProbs2(j,t) .* ProdTerm;
        end
        for i=NonEmptyInds
            ShortGammaVec(i) = sum(TempGammaProbs(IndsTabAB{i},t));
        end
        [Dummy delta2(t)] = max(ShortGammaVec);
    end
    VitStruct.AB_copy2 = delta2-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Extract all relevant parts
    VitStruct.alpha_copy = mod(VitStruct.joint_copy2, x_dim);
    VitStruct.beta_copy = (VitStruct.joint_copy2-VitStruct.alpha_copy)/x_dim;
    VitStruct.alpha_genotype = bitget(VitStruct.joint_genotype,1);
    VitStruct.beta_genotype = bitget(VitStruct.joint_genotype,2);
    % This is new and may work ...
    VitStruct.A_copy = mod(VitStruct.AB_copy2, x_dim_eq);
    VitStruct.B_copy = floor(VitStruct.AB_copy2/ x_dim_eq);
    VitStruct.A_copy = VitStruct.A_copy';
    VitStruct.B_copy = VitStruct.B_copy';

    VitStruct.total_copy = VitStruct.alpha_copy' + VitStruct.beta_copy';
    VitStruct.alpha_copy = VitStruct.alpha_copy';
    VitStruct.beta_copy = VitStruct.beta_copy';
%%    VitStruct.joint_copy = reshape(VitStruct.joint_copy(:,[1],:),x_dim*x_dim,n)';

    f{1} = 'joint_copy2';
    f{2} = 'AB_copy2';
%%    f{3} = 'joint_copy_couples';
%%    f{4} = 'AB_copy_couples';
    VitStruct.joint_copy = VitStruct.joint_copy2';
    VitStruct.AB_copy = VitStruct.AB_copy2';
    %%VitStruct.joint_copy = max(VitStruct.joint_copy')';
    VitStruct = rmfield(VitStruct, f);

end % if do_couples ...

