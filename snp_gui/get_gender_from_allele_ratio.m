%function [gender_cell, female_inds] = get_gender_from_allele_ratio(allele_ratio_vec_in_chr_x_mat, chip)
function [gender_cell, female_inds] = get_gender_from_allele_ratio(allele_ratio_vec_in_chr_x_mat, chip)

num_snps_chr_x = size(allele_ratio_vec_in_chr_x_mat, 1);
num_samples = size(allele_ratio_vec_in_chr_x_mat,2);
gender_cell = cell(num_samples, 1);
gender_cell(:) = {'M'};

temp_mat = zeros(size(allele_ratio_vec_in_chr_x_mat));
if(strcmp(lower(chip), 'xba'))
    fem_ratio_min = 0.65;
    fem_ratio_max = 2-fem_ratio_min;
    fem_frac_thresh = 0.0703;
elseif(strcmp(lower(chip), 'hind'))
    fem_ratio_min = 0.75;
    fem_ratio_max = 2-fem_ratio_min;
    fem_frac_thresh = 0.03;
end


temp_mat(find(allele_ratio_vec_in_chr_x_mat> fem_ratio_min & allele_ratio_vec_in_chr_x_mat < fem_ratio_max)) = 1;
sum_vec = sum(temp_mat);

female_inds = find(sum_vec./num_snps_chr_x > fem_frac_thresh);

if(size(female_inds,1)~= 1) female_inds = female_inds'; end
gender_cell(female_inds) = {'F'};
