% The function gives penalty to the model size. It is used for model size
% determination in the function HMM_K_by_penalized_dist_EM
function ret = penalized_dist_HMM(dist_params_vec, empirical_data)

[trans_prob_mat, dist_params] = HMM_dist_gaussian_params_vec_to_params(dist_params_vec);
num_emp_obs = length(empirical_data);

stationary_dist = get_stationary_dist(trans_prob_mat);
num_rounds = 50;
dKS = dKS_HMM_gaussian2_scan(trans_prob_mat, dist_params, empirical_data, num_rounds);

Cn = 0.01*(num_emp_obs^(-0.5))*log(num_emp_obs);

ret = dKS - Cn*sum(log(stationary_dist));




