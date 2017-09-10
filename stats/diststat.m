% Compute mean and var of a general distribution 
%
% Input:
% dist_str - string representing distribution
% num_grid_points - (optional) determine resolution of computation
%
% Output:
% mu - mean of distribution 
% V - variance of distribution 
%
function [mu V] = diststat(dist_str, num_grid_points)

if(~exist('num_grid_points', 'var')) % set default resolution
    num_grid_points = 1000;
end
res = 1/num_grid_points;
min_F = icdf(dist_str, res/1000, 0, 1); % , min(alpha)+1); % minimum possible value of distribution
max_F = icdf(dist_str, 1-res/1000, 0, 1); % , max(alpha)+1); % maximum possible value of distribution

x_vec = min_F:res:max_F;
point_mean_vec = res * pdf(dist_str, x_vec, 0, 1) .* x_vec;
point_mean_sqr_vec = point_mean_vec .* x_vec;
mu = sum(point_mean_vec); 
V = sum(point_mean_sqr_vec) - mu^2; 


