% Give hypergeometric-like p-value for pairwise intersections of many sets
%
% The function takes a number of sets and their pairwise intersections, 
% and compute the probability that the total size of all pairwise
% intersections will be larger, provided that the sets were chosen at random
% This is a generalization of the hypergeometric statistic for two sets. 
% For many sets we don't know how to compute the statistics' distribution 
% analytically. We compute the first two moments and rely on the gaussian
% approximation. 
%
% Input: 
% N - universe size. Could be either a scalar (all pairs have the same universe) 
% or a K*K matrix (different N_{ij} for different pairs. Not supported yet)
% sets_size_vec - vector of set sizes
% sets_intersection_mat - matrix of sets intersection sizes
%
% Output: 
% over_pval - prob. to get at random a larger intersection
% under_pval - prob. to get at random a smaller intersection
% z_score - the Z score which is supposed to be Normally distributed 
% mu_s - mean of s score
% sigma_s - st.d. of s score 
%
function [over_pval, under_pval, z_score, mu_s, sigma_s] = ...
    hypergeometric_for_many_sets(N, sets_size_vec, sets_intersection_mat)

K = length(sets_size_vec); % number of different sets

s = ( sum(sum(sets_intersection_mat)) - sum(diag(sets_intersection_mat)) ) / 2; % get total intersection size 

% Now compute the null model's two moments
t_over_n = sets_size_vec ./ N; 
% % % t_over_n_1 = (sets_size_vec - 1) ./ (N-1); 

t_ij_over_n_mat = t_over_n' * t_over_n;  
t_ij_over_n_mat = t_ij_over_n_mat - diag(diag(t_ij_over_n_mat));
% % % t_i_sqr_j_over_n_mat = (t_over_n').^2 * t_over_n; 
% % % t_i_sqr_j_over_n_mat = t_i_sqr_j_over_n_mat - diag(diag(t_i_sqr_j_over_n_mat));
% % % t_ij_over_n_mat_1 = t_over_n_1' * t_over_n_1;  
% % % t_ij_over_n_mat_1 = t_ij_over_n_mat_1 - diag(diag(t_ij_over_n_mat_1));
% % % t_i_i_1_j_over_n_mat = (t_over_n' .* t_over_n_1') * t_over_n; 
% % % t_i_i_1_j_over_n_mat = t_i_i_1_j_over_n_mat - diag(diag(t_i_i_1_j_over_n_mat));

t_ij_minus_n_over_n_mat = (t_over_n'-1) * (t_over_n - 1);  
t_ij_minus_n_over_n_mat = t_ij_minus_n_over_n_mat - diag(diag(t_ij_minus_n_over_n_mat));


mu_s = N * ( 0.5 * sum(t_ij_over_n_mat(:))  ); % mean of s

% % % Old way: decompose sigma to 4 
% % % sigma = zeros(4,1); 
% % % sigma(1) = mu_s - N * ( 0.5 * sum(t_ij_over_n_mat(:).^2)  ); % the var part 
% % % sigma(2) = 0.5 * N*(N-1) * sum(t_ij_over_n_mat(:) .* t_ij_over_n_mat_1(:)) - ...
% % %     0.5 * N*(N-1) * sum(t_ij_over_n_mat(:).^2);  % the cov part 1

% New way: just compute the simplified expression 
alt_sigma_1_2 = sum( t_ij_over_n_mat(:) .* t_ij_minus_n_over_n_mat(:))  * N^2 / (2*(N-1));  % a simplified expression   

% % % alt_sigma_3_4 = 0;
% % % for j=1:K % loop on the third index - but it must be different than
% the previous two !!! this loop should give zero so it's avoided
% % %     tmp_mat = t_ij_over_n_mat; tmp_mat(j,:) = 0; tmp_mat(:,j) = 0;
% % %     tmp_mat2 = t_i_sqr_j_over_n_mat; tmp_mat2(j,:) = 0; tmp_mat2(:,j) = 0;
% % %     sigma(3) = sigma(3) + N * sum(tmp_mat(:)) * t_over_n(j) - ...
% % %         N * sum(tmp_mat2(:)) * t_over_n(j); % same row overlap
% % %     tmp_mat = t_i_i_1_j_over_n_mat; tmp_mat(j,:) = 0; tmp_mat(:,j) = 0;
% % %     sigma(4) = sigma(4) + 0.5 * N*(N-1) * sum(tmp_mat(:)) * t_over_n(j) - ...
% % %         0.5 * N*(N-1) * sum(tmp_mat2(:)) * t_over_n(j); % different rows overlap
% % % 
% % %     alt_sigma_3_4 = alt_sigma_3_4 + sum(tmp_mat(:)) * t_over_n(j) * N / 2 - ...
% % %     sum(tmp_mat2(:)) * t_over_n(j) * N / 2;
% % % end


%sigma(2) = 2*sigma(2); sigma(4) = 2*sigma(4); % apply heuristic correction (needs to be changed!!!)
sigma_s = sqrt(alt_sigma_1_2); %%% sigma_s = sqrt(sum(sigma)); % take sqrt 

z_score = (s - mu_s) / sigma_s;
under_pval = normcdf( z_score ) ;% get the Zscore 
over_pval = 1-under_pval; % we ignore ties

