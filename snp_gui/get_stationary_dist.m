% This function computes the stationary distribution of 
% a Markov chain - done by taking the eigenvector 
% with eigenvalue=1 of the transposed matrix
function stationary_dist = get_stationary_dist(trans_prob_mat)

[V,D] = eig(trans_prob_mat'); % should take the eigen vector of e.value = 1;
epsilon = 10^(-5);

[min_diff min_diff_ind] = min((diag(D)-1).^2)

%% stationary_dist = V(:,find(diag(D) > 1-epsilon));
stationary_dist = V(:,min_diff_ind);
stationary_dist = stationary_dist/sum(stationary_dist);
