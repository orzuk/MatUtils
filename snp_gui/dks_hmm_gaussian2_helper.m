%function ret = dKS_HMM_gaussian2_helper(observed_data_2, trans_prob_mat,
%dist_params, empirical_dist_data)
% returns the minus of the absolute value of the difference between the two distributions on observed_data_2
function ret = dKS_HMM_gaussian2_helper(observed_data_2, trans_prob_mat, dist_params, empirical_dist_data)

ret1 = get_HMM_dist_gaussian(trans_prob_mat, dist_params, observed_data_2);
ret2 = get_empirical_dist2(empirical_dist_data, observed_data_2);

ret = -abs(ret1-ret2);
