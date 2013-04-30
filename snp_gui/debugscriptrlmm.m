%%ALL=load('F:\Or\hapmap\CEU_test\AllSamplesMat_xba.mat');
%% ALL=load('E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\CEU\AllSamplesMat_xba.mat');


load('F:\Or\SNP_GUI\database\RLMM_ALL_POPS_xba.mat');
CEU1 = load('F:\Or\SNP_GUI\database\RLMM_CEU_xba.mat');

% load('..\database\RLMM_CEU_xba.mat');
for i=1:3
    INIT_S{i} = zeros(2);
    INIT_S{i}(1,1) = RLMM.SigmaMuVec((i-1)*3+1);
    INIT_S{i}(2,2) = RLMM.SigmaMuVec((i-1)*3+2);
    INIT_S{i}(1,2) = RLMM.SigmaMuVec((i-1)*3+3); INIT_S{i}(2,1) = INIT_S{i}(1,2);
end


[temp, errors_in_rlmm] = intersect(III(Errors), I);
Errors = Errors(errors_in_rlmm);
snp_inds = [III(Errors)' setdiff(III([100:130]), III(Errors)') III(631)];
affy_snp_inds = [JJJ(Errors)' setdiff(JJJ([100:130]), JJJ(Errors)') JJJ(631)];
num_snps = length(snp_inds);

%y=III(Errors_AB(k)); %yy=J(find(I==y));
[temp_ind, temp, J_ind] = intersect_order_by_first_gr(snp_inds, I);
yy=J(J_ind);
% yy=14483; y=14573; 
k=1;
OurCall = LearnedGenotypesMat(snp_inds(k))
AffyCall = T.genotype_vec_xba(affy_snp_inds(k))
CopyMatA_all_pops = [];
CopyMatB_all_pops = [];
color_vec = 'mcg'; GenotypeMat = []; copy_num_vec = []; allele_ratio_vec = [];
hapmap_dir = 'F:\Or\hapmap\Genotype_data'; hapmap_version = genome_assembly_to_hapmap_version();
for hapmap_population = 1:3
    ALL=load(fullfile('F:\Or\hapmap\', [pop_str_vec{hapmap_population} '_test'], 'AllSamplesMat_xba.mat'));
    [ALL.NormalizedSNPsCopyMatA ALL.NormalizedSNPsCopyMatB] = ...
        replace_AB_according_to_strand(ALL.NormalizedSNPsCopyMatA, ALL.NormalizedSNPsCopyMatB, ALL.snp_ids, ALL.chip);
    [ALL_S ALL_J ALL_I ] = intersect_order_by_first_gr(RLMM.snp_ids(yy), ALL.snp_ids);

    CopyMatA_all_pops = [CopyMatA_all_pops ALL.NormalizedSNPsCopyMatA(ALL_I,:)];
    CopyMatB_all_pops = [CopyMatB_all_pops ALL.NormalizedSNPsCopyMatB(ALL_I,:)];
    ALL.snp_ids = ALL.snp_ids(ALL_I);
    
    R=load( fullfile(hapmap_dir, pop_str_vec{hapmap_population}, HAPMAP_VERSION, ...
        [chip_type, '_genotypes_chr_' pop_str_vec{hapmap_population} '_' HAPMAP_VERSION '_nr_fwd.mat']), ...
        'SnpsData', 'SnpsIDs', 'SnpsBadCalls' );

    GenotypeMat = [GenotypeMat GetSNPsGenotypeMatFromPackedR(R, RLMM.snp_ids(yy)) ];

    [tmp_copy_num_vec, tmp_allele_ratio_vec] = CopyMatsToRatio(ALL.NormalizedSNPsCopyMatA(ALL_I,:), ALL.NormalizedSNPsCopyMatB(ALL_I,:));
    
    copy_num_vec = [copy_num_vec tmp_copy_num_vec];
    allele_ratio_vec = [allele_ratio_vec tmp_allele_ratio_vec];
    
end
% choose a specific snp
[AA, AB, BB, NoCall, call_cell] = genotype_call_into_num();
figure; hold on; k= 14;

PlotAlleleRatios(copy_num_vec(k,:), allele_ratio_vec(k,:), GenotypeMat(k,:), '', 0); % No new figure
DrawRLMMOneSNP(RLMM, yy(k),1, 0); DrawRLMMOneSNP(CEU1.RLMM, yy(k),1, 0, 'mck'); 
DrawRLMMOneSNP(RLMM, -1,1, 0, '---');

plot(CopyMatA(snp_inds(k)), CopyMatB(snp_inds(k)), 'r*');
plot(CopyMatA(snp_inds(k)), CopyMatB(snp_inds(k)), 'r+');
plot(CopyMatA(snp_inds(k)), CopyMatB(snp_inds(k)), 'ro');

k
snp_inds(k)
OurCall = LearnedGenotypesMat(snp_inds(k))
AffyCall = T.genotype_vec_xba(affy_snp_inds(k))
RLMM.snp_ids{yy(k)}
title([RLMM.snp_ids{yy(k)} ' Our Call: ' call_cell{OurCall} ' AffyCall: ' call_cell{AffyCall}]);

GoodGenotypeMat = GetSNPsGenotypeMatFromPackedR(R, RLMM.snp_ids(1:10));

figure; PlotAlleleRatios(copy_num_vec(k,1:90), 1./allele_ratio_vec(k,1:90), GenotypeMat(k,1:90), '', 1); %  new figure



MixtureOfGaussiansDraw2dGaussians(vec_into_mat(RLMM.MuMuVec, 2)', INIT_S, labels_vec, legends_vec)


% Test the new Gaussian matrix method .. 
[ALL_S ALL_I ALL_J] = intersect(ALL_MAT.snp_ids, RLMM.snp_ids(yy));
X = ALL_MAT.data_A(:,ALL_I);  Y = ALL_MAT.data_B(:,ALL_I); 
% X = ALL.NormalizedSNPsCopyMatB(ALL_I,:);  Y = ALL.NormalizedSNPsCopyMatA(ALL_I,:); 
X = X(find(ALL_MAT.genotypes_AB(:,ALL_I) == AA)); Y = Y(find(ALL_MAT.genotypes_AB(:,ALL_I) == AA));


%find(Y < 0.6)); Y =Y(find(Y < 0.6));
figure; 
plot(X, Y, 'c*'); hold on;

C(1,1) = mad(X,1)*mad(X,1) * MAD_CONST_SQR;
C(2,2) = mad(Y,1)*mad(Y,1) * MAD_CONST_SQR;
C(1,2) = median( (X-median(X)) .* (Y-median(Y)) ) * MAD_CONST_SQR;
C(2,1) = C(1,2);

[U,V] = eig(C)
rotated_data = [X-median(X) Y-median(Y)] * U;
W(1,1) = mad(rotated_data(:,1),1)*mad(rotated_data(:,1),1) * MAD_CONST_SQR;
W(2,2) = mad(rotated_data(:,2),1)*mad(rotated_data(:,2),1) * MAD_CONST_SQR;
C2 = U*W*inv(U);


figure; hold on;
plot(X-median(X), Y-median(Y), 'r.');
plot(rotated_data(:,1), rotated_data(:,2), 'ko'); 


CC{1} = RLMM_AA; CC{2} = C2;
legends_vec = {'C', 'C2'};
M = zeros(2);

%M = [median(X) median(Y);median(X) median(Y)]; 
MixtureOfGaussiansDraw2dGaussians(M, CC, {'A Int.', 'B Int.'}, legends_vec);



RLMM_AA = [RLMM.SigmaMats.AA(yy,1) RLMM.SigmaMats.AA(yy,3); RLMM.SigmaMats.AA(yy,3) RLMM.SigmaMats.AA(yy,2)];
RLMM_BB = [RLMM.SigmaMats.BB(yy,1) RLMM.SigmaMats.BB(yy,3); RLMM.SigmaMats.BB(yy,3) RLMM.SigmaMats.BB(yy,2)];

det(RLMM_AA)
det(RLMM_BB)





