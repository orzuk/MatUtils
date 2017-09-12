
load('E:\Libi\tools\SNP_tool\data\Leukemia\hmm_out\HD_62_3B_d_hind_hmm_out.mat');
load('E:\Libi\tools\SNP_tool\data\Leukemia\HD_62_3B_d_hind.mat');
database = load('E:\Libi\tools\SNP_tool\database\Hind_annot_data_hg17.mat');

[C, IA, IB] = intersect(database.snp_ids, snp_ids);
chr_vec = database.chr_vec(IA);
copy_num_vec = copy_num_vec(IB);
allele_ratio_vec = allele_ratio_vec(IB);

copy_A = copy_num_vec ./ (allele_ratio_vec+1);
copy_B = copy_A .*  allele_ratio_vec;

figure;
chr10_ind = find(chr_vec==10);
plot(copy_A(chr10_ind), copy_B(chr10_ind), '.');

hold on;
chr3_ind = find(chr_vec==3);
plot(copy_A(chr3_ind), copy_B(chr3_ind), 'g.');
chr2_ind = find(chr_vec==2);
plot(copy_A(chr2_ind), copy_B(chr2_ind), 'r.');



IND_VEC = [repmat([0:4]', 5, 1) mat_into_vec(repmat([0:4]', 1, 5)')];
MEW_VEC = [repmat(HMMOutStruct.ModelParams.MEW, 5, 1) mat_into_vec(repmat(HMMOutStruct.ModelParams.MEW, 1, 5)')];
SIGMA_VEC = [repmat(HMMOutStruct.ModelParams.SIGMA, 5, 1) mat_into_vec(repmat(HMMOutStruct.ModelParams.SIGMA, 1, 5)')];

BadInds = union(find(sum(IND_VEC,2)>4), find(IND_VEC(:,1).*IND_VEC(:,2)==3));
GoodInds = setdiff(1:length(IND_VEC), BadInds);

MEW_VEC = MEW_VEC(GoodInds,:);
SIGMA_VEC = SIGMA_VEC(GoodInds,:);


labels_vec{1} ='xxx'
labels_vec{2} ='yyy'
for i=1:length(SIGMA_VEC)
    legends_vec{i} = 'x'
    SIGMA_MATS{i} = diag(SIGMA_VEC(i,:));
end
MixtureOfGaussiansDraw2dGaussians(MEW_VEC, SIGMA_MATS, labels_vec, legends_vec)

axis([-0.5 4 -0.5 4]) 
