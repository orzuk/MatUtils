%function [hmm_out, good_ind] = remove_outliers_hmm_out(hmm_out, outlier_thresh)
function [hmm_out, good_ind] = remove_outliers_hmm_out(hmm_out, outlier_thresh)

good_ind = find(hmm_out.array_outlier_vec<outlier_thresh);

hmm_out.genotype_mat = hmm_out.genotype_mat(:, good_ind);
hmm_out.average_copy_num_mat = hmm_out.average_copy_num_mat(:, good_ind);
hmm_out.copy_num_mat = hmm_out.copy_num_mat(:, good_ind);
hmm_out.zero_one_mat_del = hmm_out.zero_one_mat_del(:, good_ind);
hmm_out.zero_one_mat_amp = hmm_out.zero_one_mat_amp(:, good_ind);
hmm_out.array_outlier_vec = hmm_out.array_outlier_vec(good_ind);
