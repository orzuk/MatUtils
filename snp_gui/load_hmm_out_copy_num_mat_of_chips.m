function [err, genotype_mat, data_snp_ids, chr_num_snps_vec, average_copy_num_mat, copy_num_mat] = ...
    load_hmm_out_copy_num_mat_of_chips(in_sample_names, user_dir, chip_type, chip_snp_ids_ordered)

err='';

if strcmp(chip_type,'Hind_and_Xba')
    % load genotype data
    [genotype_mat, chr_num_snps_vec, data_snp_ids] = ...
        load_hmm_genotype_mat(in_sample_names, user_dir, 'Hind');
    [genotype_mat2, chr_num_snps_vec2, data_snp_ids2] = ...
        load_hmm_genotype_mat(in_sample_names, user_dir, 'Xba');
    chr_num_snps_vec=chr_num_snps_vec + chr_num_snps_vec2;
    genotype_mat=[genotype_mat; genotype_mat2];
    clear genotype_mat2;
    data_snp_ids=[data_snp_ids; data_snp_ids2];
    [dum idx]= ismember(chip_snp_ids_ordered, data_snp_ids);
    idx=idx(dum);
    genotype_mat=genotype_mat(idx,:);
    data_snp_ids=data_snp_ids(idx);
    if(nargout>4) % load average_copy_num_mat
        average_copy_num_mat = load_hmm_avg_copy_mat(in_sample_names, user_dir, 'Hind');
        average_copy_num_mat2 = load_hmm_avg_copy_mat(in_sample_names, user_dir, 'Xba');
        average_copy_num_mat=[average_copy_num_mat; average_copy_num_mat2];
        clear average_copy_num_mat2;
        average_copy_num_mat=average_copy_num_mat(idx,:);
    end
    if(nargout>5) % load HMM copy_num_mat
        copy_num_mat = load_hmm_copy_num_mat(in_sample_names, user_dir, 'Hind');
        copy_num_mat2 = load_hmm_copy_num_mat(in_sample_names, user_dir, 'Xba');
        copy_num_mat=[copy_num_mat; copy_num_mat2];
        clear copy_num_mat2;
        copy_num_mat=copy_num_mat(idx,:);
    end

else
    % load genotype data
    [genotype_mat, chr_num_snps_vec, data_snp_ids] = ...
        load_hmm_genotype_mat(in_sample_names, user_dir, chip_type);
    if(nargout>4) % load average_copy_num_mat
        average_copy_num_mat = load_hmm_avg_copy_mat(in_sample_names, user_dir, chip_type);
    end
    if(nargout>5) % load HMM copy_num_mat
        copy_num_mat = load_hmm_copy_num_mat(in_sample_names, user_dir, chip_type);
    end
end
