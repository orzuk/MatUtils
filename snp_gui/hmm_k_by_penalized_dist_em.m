%function [K_min, pen_dist_vec, params_cell] = HMM_K_by_penalized_dist_EM(empirical_data, K_bound)
% Gets as Input:
% 1) empirical_data - a vector of observed data
% 2) K_bound - a vector with values of K to check
%  
% output:
% 1) K_min - the value of K (out of K_bound) which gets the minimal score (here it is the best score)
% 2) pen_dist_vec - vector of scores, for each value of K
% 3) params_cell - a cell of size (length(K_bound),3). in each row it contains 3 cells:
% a) transition matrix, b) mean and std for each state (in each line mean is the first, std the second),
% c) stationary distribution

function [K_min, pen_dist_vec, params_cell, log_scores] = HMM_K_by_penalized_dist_EM(empirical_data, K_bound, num_EM_iters, num_EM_starting_points, EM_tolerance)

num_K = length(K_bound);
pen_dist_vec = zeros(1, num_K);
pen_dist_vec(:) = inf;
K_vec = K_bound;
params_cell = cell(num_K,3);

for i = 1:num_K
    Trying_Size = i
    HMM_chrom_loc = ones(size(empirical_data));
    use_locations = 0;
    HMM_x_dim = K_vec(i);
    HMM_y_dim = 1;
    place_flag = 0;    do_fold_change = 0;      mean_vec_rep = 0;      std_vec_rep = 0;
%    num_EM_iters = 100;
%    num_EM_starting_points = 10;
%    EM_tolerance = 10^(-6);
%    path(path, 'E:\Libi\tools\HMM');
    [stationary_dist trans_prob_mat mix_mat mean_mat std_mat log_score] = ...
        TrainHMMFromDataEMMatlab(empirical_data, HMM_chrom_loc, use_locations, HMM_x_dim, HMM_y_dim, ...
        place_flag, do_fold_change, mean_vec_rep, std_vec_rep, ...
        num_EM_iters, num_EM_starting_points, EM_tolerance);

    log_scores(i) = log_score;   % Save now also the score !!! 
    
    % assign values to dist_params
    dist_params = zeros(HMM_x_dim, 2);
    dist_params(:,1) = mean_mat;
    dist_params(:,2) = std_mat;
    % calc penalized score
    params_vec = HMM_dist_gaussian_params_to_vec(trans_prob_mat, dist_params);
    val = penalized_dist_HMM(params_vec, empirical_data);
    
    stationary_dist = get_stationary_dist(trans_prob_mat);
    pen_dist_vec(i) = val;
    params_cell{i,1} = trans_prob_mat;
    params_cell{i,2} = dist_params;
    params_cell{i,3} = stationary_dist;
end

[min_pen_dist min_ind] =  min(pen_dist_vec);
K_min = K_vec(min_ind)



