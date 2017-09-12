function [err, raw_copy_num_mat, data_snp_ids] = load_normalized_copy_num_mat_of_chips(sample_names, user_dir, chip_type, chip_snp_ids_ordered)

if strcmp(chip_type,'Hind_and_Xba')
    [err, raw_copy_num_mat, data_snp_ids] = load_normalized_copy_num_mat_of_chips(sample_names, user_dir, 'Hind');
    [err, raw_copy_num_mat2, data_snp_ids2] = load_normalized_copy_num_mat_of_chips(sample_names, user_dir, 'Xba');
    
    raw_copy_num_mat=[raw_copy_num_mat; raw_copy_num_mat2];
    data_snp_ids=[data_snp_ids; data_snp_ids2];
    
    [dum idx]= ismember(chip_snp_ids_ordered, data_snp_ids);
    idx=idx(dum);
    raw_copy_num_mat=raw_copy_num_mat(idx,:);
    data_snp_ids=data_snp_ids(idx);
    
else
    [err, raw_copy_num_mat, data_snp_ids] = load_normalized_copy_num_mat(sample_names, user_dir, chip_type);
end