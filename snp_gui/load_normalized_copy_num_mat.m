%function [err, raw_copy_num_mat, data_snp_ids, array_outlier_vec] = load_normalized_copy_num_mat(sample_names, user_dir, chip_type)
% loads normalized raw copy number (before the HMM ran on it)
function [err, raw_copy_num_mat, data_snp_ids, array_outlier_vec] = load_normalized_copy_num_mat(sample_names, user_dir, chip_type)

err = '';
num_samples = length(sample_names);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load normalized copy nums
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm_f = fullfile(user_dir,['AllSamplesMat_' chip_type '.mat']);

data_struct = load(norm_f);
data_sample_names = data_struct.SampleNames;
raw_copy_num_mat = data_struct.NormalizedSNPsCopyMatB+data_struct.NormalizedSNPsCopyMatA;
data_snp_ids = data_struct.snp_ids;
array_outlier_vec = data_struct.array_outlier_vec;

% leave only wanted samples
[C, IA ,IB] = intersect_order_by_first_gr(sample_names, data_sample_names);
raw_copy_num_mat = raw_copy_num_mat(:, IB);
array_outlier_vec = array_outlier_vec(IB);

%remove blank SNPs 21/03/2007
idx=find(strcmp(deblank(data_snp_ids),''));
raw_copy_num_mat(idx,:)=[];
data_snp_ids(idx)=[];
%%%%%%%%%%%%%%%%%
