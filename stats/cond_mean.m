% Compute conditional mean of a distribution given that value is > alpha
%
% Input:
% payoff_dist - string representing distribution
% alpha - value above which we take conditional mean
% num_grid_points - (optional) determine resolution of computation
%
% Output:
% cond_mean_vec - conditional expectation given we're at least alpha
%
function cond_mean_vec = cond_mean(payoff_dist, alpha, num_grid_points)

if(~exist('num_grid_points', 'var')) % set default resolution
    num_grid_points = 1000;
end
res = 1/num_grid_points;
alpha_min = min(alpha); % take lowest when computing for a vector
max_F = max(icdf(payoff_dist, 1-res/1000, 0, 1), min(max(alpha), icdf(payoff_dist, 1, 0, 1))); % maximum possible value of distribution

if(max_F == alpha) % not good - we want to take maximum above alpha 
    if(isinf(icdf(payoff_dist, 1, 0, 1)))
        max_F = alpha+0.00001;
    else
        max_F = max(alpha + 0.0000000000000000001, icdf(payoff_dist, 1, 0, 1)); % Add just  alittle bit for a bounded 
    end
end    

res = (max_F-alpha_min)/num_grid_points;
x_vec = alpha_min:res:max_F;
point_mean_vec = res * pdf(payoff_dist, alpha_min:res:max_F, 0, 1) .* x_vec;

cond_mean_vec = cumsum(point_mean_vec(end:-1:1)); % copmute cumsum from right to left
cond_mean_vec = cond_mean_vec(end:-1:1); % flip back
epsilon = 0.000000000000000000000001;
[~, ~, intersect_inds1, intersect_inds2] = ...
    intervals_intersect(alpha, alpha, x_vec, [x_vec(2:end) x_vec(end)+res], [], 1); % assume all is sorted
[intersect_inds1_unique u_inds] = unique(intersect_inds1);
%intersect_inds2(u_inds)
cond_mean_vec = cond_mean_vec(intersect_inds2(u_inds)) ./ (1-cdf(payoff_dist, alpha, 0, 1));

