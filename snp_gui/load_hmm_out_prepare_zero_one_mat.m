%function [err hmm_out] = load_hmm_out_prepare_zero_one_mat(sample_names, user_dir, chip_type, thresh_del, thresh_amp, snp_ids, chr_vec, gender)
function [err hmm_out] = load_hmm_out_prepare_zero_one_mat(sample_names, user_dir, chip_type, thresh_del, thresh_amp, snp_ids, chr_vec, gender)

[err, hmm_out.copy_num_mat, hmm_out.average_copy_num_mat, hmm_out.chr_num_snps_vec, hmm_out.genotype_mat, ...
    hmm_out.data_snp_ids] = load_hmm_out_copy_num_mat(sample_names, user_dir, chip_type);
hmm_out.genotype_mat = geno_hmm_into_affy(hmm_out.genotype_mat, hmm_out.data_snp_ids, chip_type);

[err2, hmm_out.raw_copy_num_mat, data_snp_ids, array_outlier_vec] = load_normalized_copy_num_mat(sample_names, user_dir, chip_type);

% intersect snp ids
[C, IA, IB] = intersect_order_by_first_gr(hmm_out.data_snp_ids, data_snp_ids);

hmm_out.raw_copy_num_mat = hmm_out.raw_copy_num_mat(IB, :);
hmm_out.array_outlier_vec = array_outlier_vec;

% intersect with all snp ids in order to find chrX snps
[C, IA, IB] = intersect_order_by_first_gr(hmm_out.data_snp_ids, snp_ids);
snp_chr_vec = chr_vec(IB);
chr_x_ind = find(snp_chr_vec==23);
num_snps = length(hmm_out.data_snp_ids);
non_chr_x_ind = [1:num_snps]; non_chr_x_ind(chr_x_ind) = [];

thresh_according_to_raw_flag = 1;
%thresh_according_to_raw_flag = 0;
if(thresh_according_to_raw_flag)
    hmm_out.average_copy_num_mat = hmm_out.raw_copy_num_mat;
    hmm_out.raw_copy_num_mat = []; % clear doesn't work
    del_amp_flag = 1; % for deletions analysis
    hmm_out.zero_one_mat_del = put_one_in_putative_aberr(del_amp_flag, ...
        hmm_out.average_copy_num_mat, hmm_out.chr_num_snps_vec, ...
        hmm_out.copy_num_mat, thresh_del, thresh_amp, gender, chr_x_ind);

    del_amp_flag = 2; % for amplifications analysis
    hmm_out.zero_one_mat_amp = put_one_in_putative_aberr(del_amp_flag, ...
        hmm_out.average_copy_num_mat, hmm_out.chr_num_snps_vec, ...
        hmm_out.copy_num_mat, thresh_del, thresh_amp, gender, chr_x_ind);
else
    hmm_out.raw_copy_num_mat = []; % clear doesn't work
    male_indices = strmatch('M', gender);
    female_indices = strmatch('F', gender);
    % prepare the male matrix
    zero_one_mat_del_male_non_x = zeros(length(non_chr_x_ind), length(male_indices), 'single');
    zero_one_mat_del_male_x = zeros(length(chr_x_ind), length(male_indices), 'single');
    zero_one_mat_del_male_x(find(hmm_out.average_copy_num_mat(chr_x_ind, male_indices)<=thresh_del-1))=1;
    zero_one_mat_del_male_non_x(find(hmm_out.average_copy_num_mat(non_chr_x_ind, male_indices)<=thresh_del))=1;
    zero_one_mat_del_male = zeros(num_snps, length(male_indices), 'single');
    zero_one_mat_del_male(non_chr_x_ind,:) = zero_one_mat_del_male_non_x;
    zero_one_mat_del_male(chr_x_ind,:) = zero_one_mat_del_male_x;

    zero_one_mat_amp_male_non_x = zeros(length(non_chr_x_ind), length(male_indices), 'single');
    zero_one_mat_amp_male_x = zeros(length(chr_x_ind), length(male_indices), 'single');
    zero_one_mat_amp_male_x(find(hmm_out.average_copy_num_mat(chr_x_ind, male_indices)>=thresh_amp-1))=1;
    zero_one_mat_amp_male_non_x(find(hmm_out.average_copy_num_mat(non_chr_x_ind, male_indices)>=thresh_amp))=1;
    zero_one_mat_amp_male = zeros(num_snps, length(male_indices), 'single');
    zero_one_mat_amp_male(non_chr_x_ind,:) = zero_one_mat_amp_male_non_x;
    zero_one_mat_amp_male(chr_x_ind,:) = zero_one_mat_amp_male_x;


    % prepare the female matrix
    zero_one_mat_del_female = zeros(num_snps, length(female_indices), 'single');
    zero_one_mat_del_female(find(hmm_out.average_copy_num_mat(:, female_indices)<=thresh_del))=1;
    zero_one_mat_amp_female = zeros(num_snps, length(female_indices), 'single');
    zero_one_mat_amp_female(find(hmm_out.average_copy_num_mat(:, female_indices)>=thresh_amp))=1;
    
    % combine matrics
    hmm_out.zero_one_mat_del = zeros(size(hmm_out.copy_num_mat), 'single');
    hmm_out.zero_one_mat_del(:, female_indices) = zero_one_mat_del_female;
    hmm_out.zero_one_mat_del(:, male_indices) = zero_one_mat_del_male;
    hmm_out.zero_one_mat_amp = zeros(size(hmm_out.copy_num_mat), 'single');
    hmm_out.zero_one_mat_amp(:, female_indices) = zero_one_mat_amp_female;
    hmm_out.zero_one_mat_amp(:, male_indices) = zero_one_mat_amp_male;
end
