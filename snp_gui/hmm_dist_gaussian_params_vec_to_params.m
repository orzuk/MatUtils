function [trans_prob_mat, dist_params] = HMM_dist_gaussian_params_vec_to_params(params_vec)

vec_len = length(params_vec);
init_num_states = 3;
%find number of states
num_states = fzero(inline('x.*x +2*x -a','x','a'), init_num_states, optimset('Display','off'), vec_len);
num_states = round(num_states);

trans_prob_mat = zeros(num_states, num_states);
dist_params = zeros(num_states,2);
for i = 1:num_states
    trans_prob_mat(i,:) = params_vec(i*num_states-num_states+1 :i*num_states);
end
last_ind = num_states*num_states;
dist_params(:,1) = params_vec(last_ind+1:last_ind+num_states)';
last_ind = last_ind+num_states;
dist_params(:,2) = params_vec(last_ind+1:last_ind+num_states)';
