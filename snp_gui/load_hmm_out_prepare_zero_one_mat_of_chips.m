function hmm_out=load_hmm_out_prepare_zero_one_mat_of_chips(sample_names, user_dir, chip_type, thresh_del, thresh_amp, chip_snp_ids_ordered, chr_vec, gender)

if strcmp(chip_type,'Hind_and_Xba')
%     [err, hmm_out.copy_num_mat, hmm_out.average_copy_num_mat, hmm_out.chr_num_snps_vec, hmm_out.genotype_mat, hmm_out.data_snp_ids] = ...
%         load_hmm_out_copy_num_mat(sample_names, user_dir, 'Hind');
%     hmm_out.genotype_mat = geno_hmm_into_affy(hmm_out.genotype_mat, hmm_out.data_snp_ids, 'Hind');
 
    [err hmm_out] = load_hmm_out_prepare_zero_one_mat(sample_names, user_dir, 'Hind', thresh_del, ...
        thresh_amp, chip_snp_ids_ordered, chr_vec, gender);
    [err hmm_out2] = load_hmm_out_prepare_zero_one_mat(sample_names, user_dir, 'Xba', thresh_del, ...
        thresh_amp, chip_snp_ids_ordered, chr_vec, gender);    
    
%     [err, copy_num_mat2, average_copy_num_mat2, chr_num_snps_vec2, genotype_mat2, data_snp_ids2] = load_hmm_out_copy_num_mat(sample_names, user_dir, 'Xba');
%     genotype_mat2 = geno_hmm_into_affy(genotype_mat2, data_snp_ids2, 'Xba');

    hmm_out.average_copy_num_mat=[hmm_out.average_copy_num_mat; hmm_out2.average_copy_num_mat];
    hmm_out2.average_copy_num_mat = [];
    hmm_out.copy_num_mat = [hmm_out.copy_num_mat; hmm_out2.copy_num_mat];
    hmm_out2.copy_num_mat = [];
    hmm_out.chr_num_snps_vec=hmm_out.chr_num_snps_vec + hmm_out2.chr_num_snps_vec;
    hmm_out.genotype_mat=[hmm_out.genotype_mat; hmm_out2.genotype_mat];
    hmm_out2.genotype_mat = [];
    hmm_out.data_snp_ids=[hmm_out.data_snp_ids; hmm_out2.data_snp_ids];
    hmm_out2.data_snp_ids = [];
    hmm_out.zero_one_mat_del = [hmm_out.zero_one_mat_del; hmm_out2.zero_one_mat_del];
    hmm_out2.zero_one_mat_del = [];
    hmm_out.zero_one_mat_amp = [hmm_out.zero_one_mat_amp; hmm_out2.zero_one_mat_amp];
    hmm_out2.zero_one_mat_amp = [];

    [dum idx]= ismember(chip_snp_ids_ordered, hmm_out.data_snp_ids);
    idx=idx(dum);
    hmm_out.average_copy_num_mat=hmm_out.average_copy_num_mat(idx,:);
    hmm_out.copy_num_mat=hmm_out.copy_num_mat(idx,:);
    hmm_out.genotype_mat=hmm_out.genotype_mat(idx,:);
    hmm_out.data_snp_ids=hmm_out.data_snp_ids(idx);
    hmm_out.zero_one_mat_del=hmm_out.zero_one_mat_del(idx,:);
    hmm_out.zero_one_mat_amp=hmm_out.zero_one_mat_amp(idx,:);

else
%     [err, hmm_out.copy_num_mat, hmm_out.average_copy_num_mat, hmm_out.chr_num_snps_vec, hmm_out.genotype_mat, hmm_out.data_snp_ids] = load_hmm_out_copy_num_mat(sample_names, user_dir, chip_type);
%     hmm_out.genotype_mat = geno_hmm_into_affy(hmm_out.genotype_mat, hmm_out.data_snp_ids, chip_type);
      [err hmm_out] = load_hmm_out_prepare_zero_one_mat(sample_names, user_dir, chip_type, thresh_del, thresh_amp, ...
          chip_snp_ids_ordered, chr_vec, gender);
end
