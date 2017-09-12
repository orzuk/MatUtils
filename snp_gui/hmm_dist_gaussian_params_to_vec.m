function params_vec = HMM_dist_gaussian_params_to_vec(trans_prob_mat, dist_params)

num_states = size(trans_prob_mat,1);
vec_len = num_states*num_states+2*num_states;
params_vec = zeros(1, vec_len);

for i = 1:num_states
    params_vec(i*num_states-num_states+1 :i*num_states) = trans_prob_mat(i,:);
end
last_ind = num_states*num_states;
params_vec(last_ind+1:last_ind+num_states) = dist_params(:,1)';
last_ind = last_ind+num_states;
params_vec(last_ind+1:last_ind+num_states) = dist_params(:,2)';
