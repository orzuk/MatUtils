% function data_pm = subtract_mm_from_pm(probe_intens_mat, pa_vec, ma_vec,
% pb_vec, mb_vec, mean_mm_flag) 
function data_pm = subtract_mm_from_pm(probe_intens_mat, pa_vec, ma_vec, pb_vec, mb_vec, mean_mm_flag) 

if(nargin < 6)
    mean_mm_flag = 1;
end
num_snp = size(probe_intens_mat, 1);
if(mean_mm_flag)
%    mm_avg_mat = (probe_intens_mat(:, ma_vec)+probe_intens_mat(:, mb_vec))/2;
    mm_avg_mat = add_non_nan_ind(probe_intens_mat(:, ma_vec),probe_intens_mat(:, mb_vec)); % this function adds only the non nan indices AND calculate average
    d_mat_a = probe_intens_mat(:, pa_vec)- mm_avg_mat;
    d_mat_b = probe_intens_mat(:, pb_vec)- mm_avg_mat;
else
    d_mat_a = (probe_intens_mat(:, pa_vec)-probe_intens_mat(:, ma_vec));
    d_mat_b = (probe_intens_mat(:, pb_vec)-probe_intens_mat(:, mb_vec));
end

data_pm = [d_mat_a d_mat_b];
data_pm(find(data_pm<0)) = 0;
