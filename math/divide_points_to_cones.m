% Takes as input a set of points V, and a set of pivot points 
% pivots (they may be the same points or other 
% points), and for each point finds the 'closest' pivot, where 
% closest here means the one with the least angle. the output 
% is the index of the closest pivot for each point
%
% Input: 
% V - set of points
% pivots - set of pivot points
% 
% Output: 
% cones_angle - angles to closest pivots
% cones_idx - indices of closest pivots
% 
function [cones_angle cones_idx] = divide_points_to_cones(V, pivots)

alpha = two_vecs_to_angle(V, pivots); % get angles
[cones_angle cones_idx] = min(alpha,[],2); % find minimal angles 



