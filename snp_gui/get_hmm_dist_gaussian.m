%function ret = get_HMM_dist_gaussian(trans_prob_mat, dist_params, observed_data)
% forward algorithm for HMM with Gaussian output
function ret = get_HMM_dist_gaussian(trans_prob_mat, dist_params, observed_data)

stationary_dist = get_stationary_dist(trans_prob_mat);
num_states = size(trans_prob_mat,1);

curr_forward_vec = ones(1, num_states);
% init forward values
for i = 1:num_states
    curr_forward_vec(i) = stationary_dist(i)*normcdf(observed_data(1), dist_params(i,1), dist_params(i,2));
end

num_time_points = length(observed_data);
for i = 2:num_time_points % first time point was already calculated
    prev_forward_vec = curr_forward_vec;
    for j = 1:num_states
        curr_forward_vec(j) = prev_forward_vec*trans_prob_mat(:,j)*normcdf(observed_data(i), dist_params(j,1), dist_params(j,2));
    end
end

ret = sum(curr_forward_vec);