% Compute the difference between two sets of intervals
% (everything appearing in the first set and not in the second set)
% Algorithm: compute complement of second set, and intersect with first set.
%
% Input:
% start_pos1 - starts of first interval set
% end_pos1 - ends of first interval set
% start_pos2 - starts of second interval set
% end_pos2 - ends of second interval set
% in_c_flag - use internal c function (faster. Default)
% is_sorted - assume arrays are sorted and don't sort (default false) 
% diff_mode - do we take parts of intervals (if intervals intersect) 
%
% Output:
% diff_start_pos - starts of intersection
% diff_end_pos - ends of intersection
% diff_inds1 - indices of diff intervals in first intervals
% diff_inds2 - indices of complement intervals of second intervals
%
function [diff_start_pos diff_end_pos diff_inds1 diff_inds2] = ...
    intervals_diff(start_pos1, end_pos1, start_pos2, end_pos2, in_c_flag, is_sorted, diff_mode, varargin)

if(~exist('is_sorted', 'var') || isempty(is_sorted))
    is_sorted = 0; 
end
if(~exist('in_c_flag', 'var') || isempty(in_c_flag)) % default: run in C (faster!!) 
    in_c_flag = 1;
end
FULL_INTERVALS = 1; BREAK_INTERVALS = 0; 
if(~exist('diff_mode', 'var') || isempty(diff_mode))
    diff_mode = BREAK_INTERVALS; 
end

if(diff_mode == FULL_INTERVALS) % here remove an interval even if it touches a little bit an interval from second set
    [inter_start_pos inter_end_pos inter_inds1 inter_inds2] = ... % intersect sets
        intervals_intersect(start_pos1, end_pos1, start_pos2, end_pos2, in_c_flag, is_sorted);
    diff_inds1 = setdiff(1:length(start_pos1), inter_inds1);
    diff_start_pos = start_pos1(diff_inds1); diff_end_pos = end_pos1(diff_inds1);
    diff_inds2 = ones(length(diff_inds1),1); % meaningless
else
    [comp_start_pos2 comp_end_pos2] = ... % take complement of second set
        intervals_complement(start_pos2, end_pos2, min(start_pos1), max(end_pos1), is_sorted);
    [diff_start_pos diff_end_pos diff_inds1 diff_inds2] = ... % intersect with first set
        intervals_intersect(start_pos1, end_pos1, comp_start_pos2, comp_end_pos2, in_c_flag, is_sorted);
end

pos_len_inds = find(diff_end_pos > diff_start_pos); % Remove intervals of length 0
diff_start_pos = diff_start_pos(pos_len_inds);
diff_end_pos = diff_end_pos(pos_len_inds);
diff_inds1 = diff_inds1(pos_len_inds);
diff_inds2 = diff_inds2(pos_len_inds);


