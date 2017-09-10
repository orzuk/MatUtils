%function [median_vec, mad_vec] = median_mad(mat)
function [median_vec, mad_vec] = median_mad(mat)


AssignStatsConstants();
median_vec = median(mat);
median_mat = repmat(median_vec, size(mat,1), 1);
mad_vec = MAD_CONST_SQR * median((mat-median_mat).^2);

mad_vec = sqrt(mad_vec);

