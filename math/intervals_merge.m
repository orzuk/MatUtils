% Merge (and remove overlaps) a set of intervals
%
% Input:
% start_pos - intervals starting positions
% end_pos - intervals ending positions
%
% Output:
% merged_start_pos - merged starting positions
% merged_end_pos - merged ending positions
% merge_inds - indices of original intervals in merged intervals (if more than one, peak one of originals arbitrarily)
%
function [merged_start_pos merged_end_pos merge_inds] = intervals_merge(start_pos, end_pos)

[merged_start_pos sort_perm] = sort(start_pos); merged_end_pos = end_pos(sort_perm); merge_inds = sort_perm;

new_start_pos = merged_start_pos; new_end_pos = merged_end_pos; new_inds = merge_inds; % just allocate enough size
no_overlap_flag = 0;
while(no_overlap_flag == 0)
    no_overlap_flag = 1;
    i=1; ctr=1;
    while(i < length(merged_start_pos))
        new_start_pos(ctr) = merged_start_pos(i);
        new_inds(ctr) = merge_inds(i);
        if(merged_start_pos(i+1) <= merged_end_pos(i)) % found intersection and need to unite intervals
            new_end_pos(ctr) = max(merged_end_pos(i), merged_end_pos(i+1));
            i=i+1; no_overlap_flag = 0;
        else
            new_end_pos(ctr) = merged_end_pos(i);
        end
        i=i+1; ctr=ctr+1;
    end
    if(i == length(merged_start_pos)) % just add the last segment
        new_start_pos(ctr) = merged_start_pos(i);
        new_end_pos(ctr) = merged_end_pos(i);
        new_inds(ctr) = merge_inds(i);
        i=i+1; ctr=ctr+1;
    end
    merged_start_pos = new_start_pos(1:ctr-1); merged_end_pos = new_end_pos(1:ctr-1); % update intervals
    merge_inds = new_inds(1:ctr-1); % update inds
end


