%function array_outlier_vec = calc_array_outlier(copy_num_mat, std_outlier)
function array_outlier_vec = calc_array_outlier(copy_num_mat, std_outlier)

NumSamples = size(copy_num_mat, 2);

num_snps = size(copy_num_mat,1);
std_mat = calc_std_mat(copy_num_mat);

clear copy_num_mat;

num_points = size(std_mat, 1);
% calc median std (MAD)
[median_vec, mad_vec] = median_mad(std_mat');
median_vec = median_vec';
median_mat = repmat(single(median_vec), 1, NumSamples);
mad_vec = mad_vec';
mad_mat = repmat(single(mad_vec), 1, NumSamples);
%median_dist_mat = abs(std_mat-median_mat); %wrong - need to take only high
%std sindows
median_dist_mat = std_mat-median_mat;
mad_array_outlier_mat = zeros(num_points, NumSamples, 'single');
mad_array_outlier_mat(median_dist_mat > std_outlier*mad_mat) = 1;
mad_array_outlier_vec = sum(mad_array_outlier_mat)'/num_points;


array_outlier_vec = mad_array_outlier_vec;

