% Compute the union of two sets of intervals 
% Input:
% start_pos1 - starts of first interval set
% end_pos1 - ends of first interval set
% start_pos2 - starts of second interval set
% end_pos2 - ends of second interval set
%
% Output: 
% union_start_pos - starts of union
% union_end_pos - ends of union
%
function [union_start_pos union_end_pos] = intervals_union(start_pos1, end_pos1, start_pos2, end_pos2)

if(size(start_pos1,1) == 1)
    [union_start_pos union_end_pos] = intervals_merge([start_pos1 start_pos2], [end_pos1 end_pos2]); % row vectors
else
    [union_start_pos union_end_pos] = intervals_merge([start_pos1' start_pos2'], [end_pos1' end_pos2']);
end


