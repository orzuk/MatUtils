% function [err, copy_num_mat, average_copy_num_mat, chr_num_snps_vec, genotype_mat, data_snp_ids, raw_copy_num_mat] = ...
%     load_hmm_out_copy_num_mat_and_raw(in_sample_names, user_dir, chip_type)
function [err, copy_num_mat, average_copy_num_mat, chr_num_snps_vec, genotype_mat, data_snp_ids, raw_copy_num_mat] = ...
    load_hmm_out_copy_num_mat_and_raw(in_sample_names, user_dir, chip_type)


[err, copy_num_mat, average_copy_num_mat, chr_num_snps_vec, genotype_mat, hmm_data_snp_ids] = ...
    load_hmm_out_copy_num_mat(in_sample_names, user_dir, chip_type);

[err2, raw_copy_num_mat, data_snp_ids] = load_normalized_copy_num_mat(in_sample_names, user_dir, chip_type);

% intersect snp ids

[C, IA, IB] = intersect_order_by_first_gr(hmm_data_snp_ids, data_snp_ids);

raw_copy_num_mat = raw_copy_num_mat(IB, :);
data_snp_ids = hmm_data_snp_ids;