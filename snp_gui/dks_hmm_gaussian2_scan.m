function dKS = dKS_HMM_gaussian2_scan(trans_prob_mat, dist_params, empirical_dist_data, num_rounds)

dKS = 0;

% libi: check this
observed_data_2 = zeros(1,2);
max_emp_val = max(empirical_dist_data);
min_emp_val = min(empirical_dist_data);
num_seg = sqrt(num_rounds);
seg_len = (max_emp_val-min_emp_val)/num_seg;
for i = min_emp_val-seg_len: seg_len:max_emp_val+max_emp_val
    for j = min_emp_val-seg_len: seg_len:max_emp_val+max_emp_val
        observed_data_2 = [i j];
        ret = -dKS_HMM_gaussian2_helper(observed_data_2, trans_prob_mat, dist_params, empirical_dist_data);
        if(ret > dKS)
            dKS = ret;
        end
    end
end