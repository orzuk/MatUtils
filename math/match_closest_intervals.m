% Find for each interval of a first set the closest interval of a second set 
% 
% Input: 
% start_pos1 - starts of first interval set
% end_pos1 - ends of first interval set
% start_pos2 - starts of second interval set
% end_pos2 - ends of second interval set
% 
% Output: 
% closest_ind - inds of closest elements in second set to first set
% closest_dist - distance of closest elements in second set to first set 
%
function [closest_ind closest_dist] = match_closest_intervals(start_pos1, end_pos1, start_pos2, end_pos2)

% first part is easy - find closest endpoints
[closest_ind11 closest_dist11] = match_closest_vals(start_pos1, start_pos2);
[closest_ind12 closest_dist12] = match_closest_vals(start_pos1, end_pos2);
[closest_ind21 closest_dist21] = match_closest_vals(end_pos1, start_pos2);
[closest_ind22 closest_dist22] = match_closest_vals(end_pos1, end_pos2);

[closest_dist1 closest_ind1] = min_with_inds(abs(closest_dist11), abs(closest_dist12)); 
[closest_dist2 closest_ind2] = min_with_inds(abs(closest_dist21), abs(closest_dist22)); 
[closest_dist closest_paired_ind] = min_with_inds(closest_dist1, closest_dist2); 
closest_paired_ind =closest_paired_ind .* (2+closest_ind2) + ...
    (1-closest_paired_ind) .* closest_ind1; % generate temp. inds: 0: start-start, 1:start-stop, 2:stop-start, 3:stop-stop

closest_ind = closest_ind11; 
closest_ind(closest_paired_ind == 1) = closest_ind12(closest_paired_ind == 1);
closest_ind(closest_paired_ind == 2) = closest_ind21(closest_paired_ind == 2);
closest_ind(closest_paired_ind == 3) = closest_ind22(closest_paired_ind == 3);

closest_dist = closest_dist11; 
closest_dist(closest_paired_ind == 1) = closest_dist12(closest_paired_ind == 1);
closest_dist(closest_paired_ind == 2) = closest_dist21(closest_paired_ind == 2);
closest_dist(closest_paired_ind == 3) = closest_dist22(closest_paired_ind == 3);


closest_dist((closest_paired_ind <= 1) & (closest_dist < 0)) = 0; % identify lines which fall within intervals
closest_dist((closest_paired_ind >= 2) & (closest_dist > 0)) = 0; % identify lines which fall within intervals









