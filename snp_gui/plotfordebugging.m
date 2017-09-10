%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here is an auxillary function for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotForDebugging(CHROM_MATS, user_dir, sample_name, chip_type, ...
    SNPChipAnnotStruct, HMMParamsStruct, HMM_MODEL, VitStruct, ProbsStruct, ALL_MAT)

colorvec = 'bgrkmcxo:.';     j=1;
load([user_dir '\' sample_name '_' chip_type '.mat']); % load again for the plot
snp_id_str = ['snp_id_' lower(chip_type)];
eval(['snp_id_chip = ' snp_id_str ';']);
sample_ratio_str = ['allele_ratio_vec_' lower(chip_type)];
sample_copy_num_str = ['copy_num_vec_' lower(chip_type)];
sample_genotype_call_str = ['genotype_vec_' lower(chip_type)];
eval([sample_ratio_str '= min(' sample_ratio_str ', 9999999999);']);
[SnpsNames I J] = intersect(SNPChipAnnotStruct.snp_ids, snp_id_chip); % No relation to chromosome here
HT_SIGNS = zeros(1,length(I));
HT_SIGNS(strmatch('-', SNPChipAnnotStruct.strand(I))) = 1;

% Set the relevant string
eval(['genotype_vec_' chip_type '(J(find(HT_SIGNS))) = 4-genotype_vec_' chip_type '(J(find(HT_SIGNS)));']);
eval(['genotype_vec_' chip_type '(find(genotype_vec_' chip_type ' == 0)) = 4;']);
%%        genotype_vec(J(find(HT_SIGNS))) = 4-genotype_vec(J(find(HT_SIGNS)));
%%        genotype_vec(find(genotype_vec == 0)) = 4; % Correct for the 4->0 transition

eval(['HMM_data_B = ' sample_copy_num_str './ (' sample_ratio_str ' + 1);']);
eval(['HMM_data_A = HMM_data_B .* ' sample_ratio_str ';']);
HMM_data_B = HMM_data_B(J); HMM_data_A = HMM_data_A(J);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot whole data for debugging: (should be removed later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; subplot(2,3,1); hold on;
[LLL LLL_I LLL_J] = intersect(CHROM_MATS.loc{j}, SNPChipAnnotStruct.chr_loc_vec(I)); % Get the locations of the SNPs we still have
eval(['[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(' sample_genotype_call_str '(J(LLL_J)));']);  % Get the Genotypes
plot(CHROM_MATS.data_A{j}(LLL_I(AA_ind1)),CHROM_MATS.data_B{j}(LLL_I(AA_ind1)), '.');
plot(CHROM_MATS.data_A{j}(LLL_I(AB_ind1)),CHROM_MATS.data_B{j}(LLL_I(AB_ind1)), 'r.');
plot(CHROM_MATS.data_A{j}(LLL_I(BB_ind1)),CHROM_MATS.data_B{j}(LLL_I(BB_ind1)), 'g.');
plot(CHROM_MATS.data_A{j}(LLL_I(no_call_ind1)),CHROM_MATS.data_B{j}(LLL_I(no_call_ind1)), 'm.');
legend('AA', 'AB', 'BB', 'no call'); xlabel('A Intensity'); ylabel('B Intensity'); title(['Sample '  sample_name  ' Chrom. ' num2str(j) ' - Affy Genotypes']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%[ALL_NC ALL_AB ALL_AA ALL_BB] = calls_vec_into_ind(genotype_vec(J));
subplot(2,3,4); hold on;
tmp_ctr = 0; ALL_AA = []; ALL_AB = []; ALL_BB = []; ALL_NC = [];
for c=HMMParamsStruct.ChromosomesToRun
    [LLL LLL_I LLL_J] = intersect(CHROM_MATS.loc{c}, SNPChipAnnotStruct.chr_loc_vec(I)); % Get the locations of the SNPs we still have
    eval(['[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(' sample_genotype_call_str '(J(LLL_J)));']);  % Get the Genotypes
    ALL_AA = [ALL_AA (LLL_I(AA_ind1)+tmp_ctr)'];
    ALL_AB = [ALL_AB (LLL_I(AB_ind1)+tmp_ctr)'];
    ALL_BB = [ALL_BB (LLL_I(BB_ind1)+tmp_ctr)'];
    ALL_NC = [ALL_NC (LLL_I(no_call_ind1)+tmp_ctr)'];
    tmp_ctr = tmp_ctr + length(VitStruct{c}.alpha_genotype);
end
plot(ALL_MAT.data_A(ALL_AA),ALL_MAT.data_B(ALL_AA), '.');
plot(ALL_MAT.data_A(ALL_AB),ALL_MAT.data_B(ALL_AB), 'r.');
plot(ALL_MAT.data_A(ALL_BB),ALL_MAT.data_B(ALL_BB), 'g.');
plot(ALL_MAT.data_A(ALL_NC),ALL_MAT.data_B(ALL_NC), 'm.');
legend('AA', 'AB', 'BB', 'no call'); xlabel('A Intensity'); ylabel('B Intensity'); title('All Chromosomes - Affy Genotypes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot First Chromosome for debugging: (should be removed later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,2); hold on;
A_ind = find(VitStruct{j}.alpha_genotype==0);  B_ind = find(VitStruct{j}.alpha_genotype);
A_ind2 = find(VitStruct{j}.beta_genotype==0);  B_ind2 = find(VitStruct{j}.beta_genotype);
BB_ind = intersect(B_ind, B_ind2);
AA_ind = intersect(A_ind, A_ind2);
AB_ind = union(intersect(A_ind, B_ind2), intersect(B_ind, A_ind2));
no_call_ind = [];
plot(CHROM_MATS.data_A{j}(AA_ind),CHROM_MATS.data_B{j}(AA_ind), '.');
plot(CHROM_MATS.data_A{j}(AB_ind),CHROM_MATS.data_B{j}(AB_ind), 'r.');
plot(CHROM_MATS.data_A{j}(BB_ind),CHROM_MATS.data_B{j}(BB_ind), 'g.');
plot(CHROM_MATS.data_A{j}(no_call_ind),CHROM_MATS.data_B{j}(no_call_ind), 'm.');
legend('AA', 'AB', 'BB'); xlabel('A Intensity'); ylabel('B Intensity'); title(['Chromosome ' num2str(j) ' - HMM Genotypes']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot whole data for debugging: (should be removed later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,5); hold on;
tmp_ctr = 0; ALL_AA = []; ALL_AB = []; ALL_BB = [];
for c=HMMParamsStruct.ChromosomesToRun
    A_ind = find(VitStruct{c}.alpha_genotype==0);  B_ind = find(VitStruct{c}.alpha_genotype);
    A_ind2 = find(VitStruct{c}.beta_genotype==0);  B_ind2 = find(VitStruct{c}.beta_genotype);
    BB_ind = intersect(B_ind, B_ind2);
    AA_ind = intersect(A_ind, A_ind2);
    AB_ind = union(intersect(A_ind, B_ind2), intersect(B_ind, A_ind2));
    ALL_AA = [ALL_AA (AA_ind+tmp_ctr)'];
    ALL_AB = [ALL_AB (AB_ind+tmp_ctr)'];
    ALL_BB = [ALL_BB (BB_ind+tmp_ctr)'];
    tmp_ctr = tmp_ctr + length(VitStruct{c}.alpha_genotype);
end
plot(ALL_MAT.data_A(ALL_AA),ALL_MAT.data_B(ALL_AA), '.');
plot(ALL_MAT.data_A(ALL_AB),ALL_MAT.data_B(ALL_AB), 'r.');
plot(ALL_MAT.data_A(ALL_BB),ALL_MAT.data_B(ALL_BB), 'g.');
legend('AA', 'AB', 'BB'); xlabel('A Intensity'); ylabel('B Intensity'); title('All Chromosomes - HMM Genotypes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(2,3,3); hold on;
[MaxMarginalGeno MaxMarginalGenoIndex] = max(ProbsStruct{j}.joint_genotype');
AA_ind = find(MaxMarginalGenoIndex==1);
AB_ind = union(find(MaxMarginalGenoIndex==2), find(MaxMarginalGenoIndex==3));
BB_ind = find(MaxMarginalGenoIndex==4);
no_call_ind = [];
plot(CHROM_MATS.data_A{j}(AA_ind),CHROM_MATS.data_B{j}(AA_ind), '.');
plot(CHROM_MATS.data_A{j}(AB_ind),CHROM_MATS.data_B{j}(AB_ind), 'r.');
plot(CHROM_MATS.data_A{j}(BB_ind),CHROM_MATS.data_B{j}(BB_ind), 'g.');
plot(CHROM_MATS.data_A{j}(no_call_ind),CHROM_MATS.data_B{j}(no_call_ind), 'm.');
legend('AA', 'AB', 'BB'); xlabel('A Intensity'); ylabel('B Intensity'); title(['Chromosome ' num2str(j) ' - HMM Genotypes Marginal']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(2,3,6); hold on;
tmp_ctr = 0; ALL_AA = []; ALL_AB = []; ALL_BB = [];
for c=HMMParamsStruct.ChromosomesToRun
    [MaxMarginalGeno MaxMarginalGenoIndex] = max(ProbsStruct{c}.joint_genotype');
    AA_ind = find(MaxMarginalGenoIndex==1);
    AB_ind = union(find(MaxMarginalGenoIndex==2), find(MaxMarginalGenoIndex==3));
    BB_ind = find(MaxMarginalGenoIndex==4);
    no_call_ind = [];
    ALL_AA = [ALL_AA (AA_ind+tmp_ctr)];
    ALL_AB = [ALL_AB (AB_ind+tmp_ctr)];
    ALL_BB = [ALL_BB (BB_ind+tmp_ctr)];
    tmp_ctr = tmp_ctr + length(VitStruct{c}.alpha_genotype);
end
plot(ALL_MAT.data_A(ALL_AA),ALL_MAT.data_B(ALL_AA), '.');
plot(ALL_MAT.data_A(ALL_AB),ALL_MAT.data_B(ALL_AB), 'r.');
plot(ALL_MAT.data_A(ALL_BB),ALL_MAT.data_B(ALL_BB), 'g.');
legend('AA', 'AB', 'BB'); xlabel('A Intensity'); ylabel('B Intensity'); title('All Chromosomes - HMM Genotypes Marginal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot whole data for debugging: (should be removed later) This time with locations !!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; subplot(2,3,1); hold on;
[LLL LLL_I LLL_J] = intersect(CHROM_MATS.loc{j}, SNPChipAnnotStruct.chr_loc_vec(I)); % Get the locations of the SNPs we still have
eval(['[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(' sample_genotype_call_str '(J(LLL_J)));']);  % Get the Genotypes
plot(CHROM_MATS.loc{j}(LLL_I(AA_ind1)),CHROM_MATS.data_A{j}(LLL_I(AA_ind1)), '.');
plot(CHROM_MATS.loc{j}(LLL_I(AB_ind1)),CHROM_MATS.data_A{j}(LLL_I(AB_ind1)), 'r.');
plot(CHROM_MATS.loc{j}(LLL_I(BB_ind1)),CHROM_MATS.data_A{j}(LLL_I(BB_ind1)), 'g.');
plot(CHROM_MATS.loc{j}(LLL_I(no_call_ind1)),CHROM_MATS.data_A{j}(LLL_I(no_call_ind1)), 'm.');
legend('AA', 'AB', 'BB', 'no call'); xlabel('Chrom. Loc.'); ylabel('A Intensity'); title(['Sample ' sample_name 'Chrom. ' num2str(j) ' - Affy Genotypes']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%[ALL_NC ALL_AB ALL_AA ALL_BB] = calls_vec_into_ind(genotype_vec(J));
subplot(2,3,4); hold on;
tmp_ctr = 0; ALL_AA = []; ALL_AB = []; ALL_BB = []; ALL_NC = [];
for c=HMMParamsStruct.ChromosomesToRun
    [LLL LLL_I LLL_J] = intersect(CHROM_MATS.loc{c}, SNPChipAnnotStruct.chr_loc_vec(I)); % Get the locations of the SNPs we still have
    eval(['[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(' sample_genotype_call_str '(J(LLL_J)));']);  % Get the Genotypes
    ALL_AA = [ALL_AA (LLL_I(AA_ind1)+tmp_ctr)'];
    ALL_AB = [ALL_AB (LLL_I(AB_ind1)+tmp_ctr)'];
    ALL_BB = [ALL_BB (LLL_I(BB_ind1)+tmp_ctr)'];
    ALL_NC = [ALL_NC (LLL_I(no_call_ind1)+tmp_ctr)'];
    tmp_ctr = tmp_ctr + length(VitStruct{c}.alpha_genotype);
end
plot(ALL_MAT.loc(ALL_AA),ALL_MAT.data_A(ALL_AA), '.');
plot(ALL_MAT.loc(ALL_AB),ALL_MAT.data_A(ALL_AB), 'r.');
plot(ALL_MAT.loc(ALL_BB),ALL_MAT.data_A(ALL_BB), 'g.');
plot(ALL_MAT.loc(ALL_NC),ALL_MAT.data_A(ALL_NC), 'm.');
legend('AA', 'AB', 'BB', 'no call'); xlabel('Chrom. Loc.'); ylabel('A Intensity'); title('All Chromosomes - Affy Genotypes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot First Chromosome for debugging: (should be removed later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,2); hold on;
A_ind = find(VitStruct{j}.alpha_genotype==0);  B_ind = find(VitStruct{j}.alpha_genotype);
A_ind2 = find(VitStruct{j}.beta_genotype==0);  B_ind2 = find(VitStruct{j}.beta_genotype);
BB_ind = intersect(B_ind, B_ind2);
AA_ind = intersect(A_ind, A_ind2);
AB_ind = union(intersect(A_ind, B_ind2), intersect(B_ind, A_ind2));
no_call_ind = [];
plot(CHROM_MATS.loc{j}(AA_ind),CHROM_MATS.data_A{j}(AA_ind), '.');
plot(CHROM_MATS.loc{j}(AB_ind),CHROM_MATS.data_A{j}(AB_ind), 'r.');
plot(CHROM_MATS.loc{j}(BB_ind),CHROM_MATS.data_A{j}(BB_ind), 'g.');
plot(CHROM_MATS.loc{j}(no_call_ind),CHROM_MATS.data_A{j}(no_call_ind), 'm.');
legend('AA', 'AB', 'BB'); xlabel('Chrom. Loc.'); ylabel('A Intensity'); title(['Chromosome ' num2str(j) ' - HMM Genotypes']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot whole data for debugging: (should be removed later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,5); hold on;
tmp_ctr = 0; ALL_AA = []; ALL_AB = []; ALL_BB = [];
for c=HMMParamsStruct.ChromosomesToRun
    A_ind = find(VitStruct{c}.alpha_genotype==0);  B_ind = find(VitStruct{c}.alpha_genotype);
    A_ind2 = find(VitStruct{c}.beta_genotype==0);  B_ind2 = find(VitStruct{c}.beta_genotype);
    BB_ind = intersect(B_ind, B_ind2);
    AA_ind = intersect(A_ind, A_ind2);
    AB_ind = union(intersect(A_ind, B_ind2), intersect(B_ind, A_ind2));
    ALL_AA = [ALL_AA (AA_ind+tmp_ctr)'];
    ALL_AB = [ALL_AB (AB_ind+tmp_ctr)'];
    ALL_BB = [ALL_BB (BB_ind+tmp_ctr)'];
    tmp_ctr = tmp_ctr + length(VitStruct{c}.alpha_genotype);
end
plot(ALL_MAT.loc(ALL_AA),ALL_MAT.data_A(ALL_AA), '.');
plot(ALL_MAT.loc(ALL_AB),ALL_MAT.data_A(ALL_AB), 'r.');
plot(ALL_MAT.loc(ALL_BB),ALL_MAT.data_A(ALL_BB), 'g.');
legend('AA', 'AB', 'BB'); xlabel('Chrom. Loc.'); ylabel('A Intensity'); title('All Chromosomes - HMM Genotypes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(2,3,3); hold on;
[MaxMarginalGeno MaxMarginalGenoIndex] = max(ProbsStruct{j}.joint_genotype');
AA_ind = find(MaxMarginalGenoIndex==1);
AB_ind = union(find(MaxMarginalGenoIndex==2), find(MaxMarginalGenoIndex==3));
BB_ind = find(MaxMarginalGenoIndex==4);
no_call_ind = [];
plot(CHROM_MATS.loc{j}(AA_ind),CHROM_MATS.data_A{j}(AA_ind), '.');
plot(CHROM_MATS.loc{j}(AB_ind),CHROM_MATS.data_A{j}(AB_ind), 'r.');
plot(CHROM_MATS.loc{j}(BB_ind),CHROM_MATS.data_A{j}(BB_ind), 'g.');
plot(CHROM_MATS.loc{j}(no_call_ind),CHROM_MATS.data_A{j}(no_call_ind), 'm.');
legend('AA', 'AB', 'BB'); xlabel('Chrom. Loc.'); ylabel('A Intensity'); title(['Chromosome ' num2str(j) ' - HMM Genotypes Marginal']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(2,3,6); hold on;
tmp_ctr = 0; ALL_AA = []; ALL_AB = []; ALL_BB = [];
for c=HMMParamsStruct.ChromosomesToRun
    [MaxMarginalGeno MaxMarginalGenoIndex] = max(ProbsStruct{c}.joint_genotype');
    AA_ind = find(MaxMarginalGenoIndex==1);
    AB_ind = union(find(MaxMarginalGenoIndex==2), find(MaxMarginalGenoIndex==3));
    BB_ind = find(MaxMarginalGenoIndex==4);
    no_call_ind = [];
    ALL_AA = [ALL_AA (AA_ind+tmp_ctr)];
    ALL_AB = [ALL_AB (AB_ind+tmp_ctr)];
    ALL_BB = [ALL_BB (BB_ind+tmp_ctr)];
    tmp_ctr = tmp_ctr + length(VitStruct{c}.alpha_genotype);
end
plot(ALL_MAT.loc(ALL_AA),ALL_MAT.data_A(ALL_AA), '.');
plot(ALL_MAT.loc(ALL_AB),ALL_MAT.data_A(ALL_AB), 'r.');
plot(ALL_MAT.loc(ALL_BB),ALL_MAT.data_A(ALL_BB), 'g.');
legend('AA', 'AB', 'BB'); xlabel('Chrom. Loc.'); ylabel('A Intensity'); title('All Chromosomes - HMM Genotypes Marginal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55




% Newest : Save in the format required by: 1. The displaying program
%                                          2. Further multi-SNP analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. The displaying program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
figure; subplot(3,3,1); hold on;
Copy_ind = {}; legend_vec = {}; cc=1;
for c=0:4
    Copy_ind{c+1} = find(VitStruct{j}.total_copy == c);
    if(~isempty(Copy_ind{c+1}))
        legend_vec{cc} = num2str(c); cc=cc+1;
    end
    plot(CHROM_MATS.data_A{j}(Copy_ind{c+1}),CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('A Intensity'); ylabel('B Intensity'); title(['Sample ' sample_name ' Chrom. ' num2str(j) ' - HMM total Copy number']);
legend(legend_vec);
%    legend('0', '1', '2', '3', '4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,2); hold on;
Copy_ind = {};
for c=0:4
    Copy_ind{c+1} = find(VitStruct{j}.A_copy == c);
    plot(CHROM_MATS.data_A{j}(Copy_ind{c+1}),CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('A Intensity'); ylabel('B Intensity'); title('Chromosome one - HMM A Copy number');
legend('0', '1', '2', '3', '4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,3); hold on;
Copy_ind = {};
for c=0:4
    Copy_ind{c+1} = find(VitStruct{j}.B_copy == c);
    plot(CHROM_MATS.data_A{j}(Copy_ind{c+1}),CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('A Intensity'); ylabel('B Intensity'); title('Chromosome one - HMM B Copy number');
legend('0', '1', '2', '3', '4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,4); hold on;
Copy_ind = {};
for c=0:2
    Copy_ind{c+1} = find(VitStruct{j}.alpha_copy == c);
    plot(CHROM_MATS.data_A{j}(Copy_ind{c+1}),CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('A Intensity'); ylabel('B Intensity'); title('Chromosome one - HMM \alpha Copy number');
legend('0', '1', '2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,5); hold on;
Copy_ind = {};
for c=0:2
    Copy_ind{c+1} = find(VitStruct{j}.beta_copy == c);
    plot(CHROM_MATS.data_A{j}(Copy_ind{c+1}),CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('A Intensity'); ylabel('B Intensity'); title('Chromosome one - HMM \beta Copy number');
legend('0', '1', '2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,6); hold on;
Copy_ind = {};
for c=0:4
    Copy_ind{c+1} = find(VitStruct{j}.A_copy +  VitStruct{j}.B_copy == c);
    plot(CHROM_MATS.loc{j}(Copy_ind{c+1}),CHROM_MATS.data_A{j}(Copy_ind{c+1})+CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('A Intensity'); ylabel('B Intensity'); title('Chromosome one - HMM Total Copy number');
legend('0', '1', '2', '3', '4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,7); hold on;
Copy_ind = {};
for c=0:2
    Copy_ind{c+1} = find(VitStruct{j}.alpha_copy == c);
    plot(CHROM_MATS.loc{j}(Copy_ind{c+1}),CHROM_MATS.data_A{j}(Copy_ind{c+1})+CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('Chrom. Loc.'); ylabel('A+B Intensity'); title('Chromosome one - HMM \alpha Copy number');
legend('0', '1', '2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,8); hold on;
Copy_ind = {};
for c=0:4
    Copy_ind{c+1} = find(VitStruct{j}.A_copy == c);
    plot(CHROM_MATS.loc{j}(Copy_ind{c+1}),CHROM_MATS.data_A{j}(Copy_ind{c+1})+CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('Chrom. Loc.'); ylabel('A+B Intensity'); title('Chromosome one - HMM A Copy number');
legend('0', '1', '2','3','4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(3,3,9); hold on;
Copy_ind = {};
for c=0:4
    Copy_ind{c+1} = find(VitStruct{j}.B_copy == c);
    plot(CHROM_MATS.loc{j}(Copy_ind{c+1}),CHROM_MATS.data_A{j}(Copy_ind{c+1})+CHROM_MATS.data_B{j}(Copy_ind{c+1}), [colorvec(c+1) '.']);
end
xlabel('Chrom. Loc.'); ylabel('A+B Intensity'); title('Chromosome one - HMM B Copy number');
legend('0', '1', '2','3','4');


%% figure; hold on;
%% plot(CHROM_MATS.loc{j},VitStruct{j}.A_copy +  VitStruct{j}.B_copy);
%% plot(CHROM_MATS.loc{j},smooth(CHROM_MATS.data_A{j}+CHROM_MATS.data_B{j}, 20), 'r.');
%% xlabel('A Intensity'); ylabel('B Intensity'); title('Chromosome one - HMM Total Copy number');
%% legend('0', '1', '2', '3', '4');

% This need not be done anymore (we've already got the files saved)
%%%    all_probe_intens_txt_into_mat_files('E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\', 'E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot classification according to the HMM
figure; hold on; plot(CHROM_MATS.data_A{j},CHROM_MATS.data_B{j}, '.')
theta = [0:0.01:2*pi]'; x = [cos(theta) sin(theta)]';
jj = 1;
for a=1:5
    for b=1:6-a  % The maximal sum is 4 copy numbers (each chromosome is either 0, 1 or 2)
        plot(HMM_MODEL.MEW(a), HMM_MODEL.MEW(b), '+k');
        plot(HMM_MODEL.MEW(a), HMM_MODEL.MEW(b), '*k');
        plot(HMM_MODEL.MEW(a), HMM_MODEL.MEW(b), 'xk');
        plot(HMM_MODEL.MEW(a), HMM_MODEL.MEW(b), 'ok');
        y(1,:) = HMM_MODEL.SIGMA(a) * x(1,:) + repmat(HMM_MODEL.MEW(a),  length(x), 1)';
        y(2,:) = HMM_MODEL.SIGMA(b) * x(2,:) + repmat(HMM_MODEL.MEW(b),  length(x), 1)';
        plot(y(1,:), y(2,:), 'k.');
        M_IND(jj,1) = HMM_MODEL.MEW(a); M_2d(jj,1) = HMM_MODEL.MEW(a); M_2d(jj,2) = HMM_MODEL.MEW(b);
        SIGMA_IND(jj,1) = HMM_MODEL.SIGMA(a); SIGMA_IND(jj,2) = HMM_MODEL.SIGMA(b);
        SIGMA_2d{jj} = [HMM_MODEL.SIGMA(a) 0; 0 HMM_MODEL.SIGMA(b)];
        P_IND(jj,1) = 1; P_IND(jj,2) = 1;
        jj = jj+1;
    end
end
P_IND = (P_IND ./ sum(sum(P_IND))) .* 15;     % Normalize Prior





















