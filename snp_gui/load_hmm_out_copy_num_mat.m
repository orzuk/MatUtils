%function [raw_copy_num_mat, chr_num_snps_vec, genotype_mat, data_snp_ids] =
%load_hmm_out_copy_num_mat2(sample_names, user_dir, chip_type, data_snp_ids)
function [err, copy_num_mat, average_copy_num_mat, chr_num_snps_vec, genotype_mat, data_snp_ids] = ...
    load_hmm_out_copy_num_mat(in_sample_names, user_dir, chip_type)

err='';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load hmm output for all samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

genotype_mat = load_hmm_genotype_mat(in_sample_names, user_dir, chip_type);
average_copy_num_mat = load_hmm_avg_copy_mat(in_sample_names, user_dir, chip_type);
copy_num_mat = load_hmm_copy_num_mat(in_sample_names, user_dir, chip_type);

f_name = fullfile(user_dir, 'display', ['all_samples_hmm_out_' chip_type '.mat']);
load(f_name, 'chr_num_snps_vec', 'data_snp_ids');
