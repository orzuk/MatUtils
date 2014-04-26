% Sample a set of genotypes and phenotype prob. for cases and controls
%
% Input:
% N - number of loci
% dispose_prob_frac - ignore a certain part of the x's probability space
% f - marginal frequencies
% compute_method_flag - how to compute (sampling, analytic etc.)
% iters - number of different genotypes to sample
%
% Output:
% x_vec - vector of genotypes (binary)
% p_x_vec - vector of probabilities Pr(x = i)
% x_ind_vec - vector of indices (???)
% x_ind_mat - matrix of ???
%
function [x_vec_cases p_x_vec_cases z_vec_cases ...
    x_vec_controls p_x_vec_controls z_vec_controls] = ...
    sample_case_control(N, dispose_prob_frac, f, ...
    compute_method_flag, num_cases, num_controls, ...
    params_struct, architecture_str)

num_assigned_cases = 0;
num_assigned_controls = 0;

iters = num_cases + num_controls;
while(num_assigned_cases < num_cases) % simulate many genotypes
    [x_vec p_x_vec] = ...
        initilize_x_vec_constants(N, dispose_prob_frac, f, ...
        compute_method_flag, iters);
    z_vec = genetic_architecture(x_vec, architecture_str, params_struct);
    z_vec_binary = z_vec > rand(iters,1); % assign disease value
    cases_inds = find(z_vec_binary == 1);
    controls_inds = find(z_vec_binary == 0);
    num_new_cases = min(length(cases_inds), ...
        num_cases - num_assigned_cases);
    x_vec_cases(num_assigned_cases+1:num_assigned_cases+num_new_cases,:) = ...
        x_vec(cases_inds(1:num_new_cases),:);
    p_x_vec_cases(num_assigned_cases+1:num_assigned_cases+num_new_cases,:) = ...
        p_x_vec(cases_inds(1:num_new_cases),:);
    z_vec_cases(num_assigned_cases+1:num_assigned_cases+num_new_cases,:) = ...
        z_vec(cases_inds(1:num_new_cases),:);
    num_new_controls = min(iters-length(cases_inds), ...
        num_controls - num_assigned_controls);
    x_vec_controls(num_assigned_controls+1:num_assigned_controls+num_new_controls,:) = ...
        x_vec(controls_inds(1:num_new_controls),:);
    p_x_vec_controls(num_assigned_controls+1:num_assigned_controls+num_new_controls,:) = ...
        p_x_vec(controls_inds(1:num_new_controls),:);
    z_vec_controls(num_assigned_controls+1:num_assigned_controls+num_new_controls,:) = ...
        z_vec(controls_inds(1:num_new_controls),:);
    num_assigned_cases = num_assigned_cases + length(num_new_cases);
    num_assigned_controls = num_assigned_controls + length(num_new_controls);
end




