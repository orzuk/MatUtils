% Compute the intersection between two sets of intervals
%
% Input:
% start_pos1 - starts of first interval set
% end_pos1 - ends of first interval set
% start_pos2 - starts of second interval set
% end_pos2 - ends of second interval set
% in_c_flag - use internal c function (faster. Default)
% is_sorted - assume arrays are sorted and don't sort (default false) 
%
% Output:
% intersect_start_pos - starts of intersection
% intersect_end_pos - ends of intersection
% intersect_inds1 - indices of intersect intervals in first intervals
% intersect_inds2 - indices of intersect intervals in second intervals
%
function [intersect_start_pos intersect_end_pos intersect_inds1 intersect_inds2] = ...
    intervals_intersect(start_pos1, end_pos1, start_pos2, end_pos2, in_c_flag, is_sorted, varargin)

if(~exist('is_sorted', 'var') || isempty(is_sorted)) % default: not sorted (need to sort) 
    is_sorted = 0; 
end
if(~exist('in_c_flag', 'var') || isempty(in_c_flag)) % default: run in C (faster!!) 
    in_c_flag = 1;
end

if(~is_sorted) % always do sorting in matlab (faster)
    [start_pos1 end_pos1 merge_inds1] = intervals_merge(start_pos1, end_pos1); % merge overlapping intervals (otherwise function may not work)
    [start_pos2 end_pos2 merge_inds2] = intervals_merge(start_pos2, end_pos2); % merge overlapping intervals (otherwise function may not work)
    
    tmp_max = max(start_pos1, end_pos1); tmp_min = min(start_pos1, end_pos1); start_pos1 = tmp_min; end_pos1 = tmp_max;
    tmp_max = max(start_pos2, end_pos2); tmp_min = min(start_pos2, end_pos2); start_pos2 = tmp_min; end_pos2 = tmp_max; % make sure end > start
    [start_pos1 sort_perm1] = sort(start_pos1);
    end_pos1 = end_pos1(sort_perm1); % sort the first positions vec by start
    [start_pos2 sort_perm2] = sort(start_pos2);
    end_pos2 = end_pos2(sort_perm2); % sort the second positions vec by start
end

if(isempty(start_pos1) || isempty(start_pos2)) % can't intersect empty intervals 
    intersect_start_pos = []; intersect_end_pos = []; intersect_inds1 = []; intersect_inds2 = []; 
    return; 
end

if(in_c_flag) % This works so far only for doubles ! (result is returned in doubles too)
    %    intervals_length_to_allocate = length(start_pos1) + length(start_pos2)
    %    memory_length_before = memory
    %    memory
    
    [intersect_start_pos intersect_end_pos intersect_inds1 intersect_inds2] = ...
        intervals_intersect_c(double(start_pos1), double(end_pos1), ...
        double(start_pos2), double(end_pos2), 0*is_sorted); % for now always sort inside 
    %    memory_length_after = memory
    %    memory
    if(isempty(intersect_start_pos))
        intersect_start_pos = intersect_start_pos';
        intersect_end_pos = intersect_end_pos';
        intersect_inds1 = intersect_inds1';
        intersect_inds2 = intersect_inds2';
    end
else
    n1 = length(start_pos1); n2 = length(start_pos2);
    
 %   sort_perm1_is = sort_perm1, sort_perm2_is = sort_perm2
    prev_ctr2 = 1; inter_ctr = 1; % counters following both regions
    intersect_start_pos = zeros(1, max(n1,n2)); intersect_end_pos = zeros(1, max(n1,n2)); 
    intersect_inds1 = zeros(1, max(n1,n2)); intersect_inds2 = zeros(1, max(n1,n2));
    for ctr1 = 1:n1 % loop on first intervals start positions 
        ctr2 = find(start_pos1(ctr1) <= end_pos2(prev_ctr2:end), 1); % find the first chance for intersection
        if(~isempty(ctr2))
            inter_start = max(start_pos1(ctr1), start_pos2(prev_ctr2+ctr2-1));
            inter_end = min(end_pos1(ctr1), end_pos2(prev_ctr2+ctr2-1));
%            inter_start_is = inter_start, inter_end_is = inter_end
%            sprintf('before while: start %f, %f\n', start_pos1(ctr1), start_pos2(prev_ctr2+ctr2-1))
%            sprintf('before while: end %f, %f\n', end_pos1(ctr1), end_pos2(prev_ctr2+ctr2-1))
%            sprintf('before while inter_start %f before while inter_end %f\n', inter_start, inter_end)
            while(inter_start <= inter_end)
                intersect_start_pos(inter_ctr) = inter_start;
                intersect_end_pos(inter_ctr) = inter_end;
                intersect_inds1(inter_ctr) = ctr1;   % New: get also the indexes
                intersect_inds2(inter_ctr) = prev_ctr2+ctr2-1;
                inter_ctr = inter_ctr+1;
                if(prev_ctr2+ctr2 <= n2) % make sure we didn't reach the end of second array
                    inter_start = max(start_pos1(ctr1), start_pos2(prev_ctr2+ctr2));
                    inter_end = min(end_pos1(ctr1), end_pos2(prev_ctr2+ctr2));
                else
                    inter_start = 1; inter_end = -1;
                end
%                sprintf('inside while inter_start %f inside while inter_end %f\n', inter_start, inter_end)
                ctr2=ctr2+1;
            end
            prev_ctr2 = max(prev_ctr2, prev_ctr2+ctr2-2); % shift previous counter to the right
%            inter_start_is2 = inter_start, inter_end_is2 = inter_end
        end
    end
%    inter_ctr_is = inter_ctr
    intersect_start_pos = intersect_start_pos(1:inter_ctr-1); intersect_end_pos = intersect_end_pos(1:inter_ctr-1); % take only relevant inds
    intersect_inds1 = intersect_inds1(1:inter_ctr-1); intersect_inds2 = intersect_inds2(1:inter_ctr-1); % take only relevant inds
    
end % if in_c_flag


if(~is_sorted)
    intersect_inds1 = sort_perm1(intersect_inds1); intersect_inds2 = sort_perm2(intersect_inds2);  % transfer back to original indices
    intersect_inds1 = merge_inds1(intersect_inds1); intersect_inds2 = merge_inds2(intersect_inds2);
end

if(isrowvector(start_pos1) && isrowvector(start_pos2))
    intersect_start_pos = vec2row(intersect_start_pos);
    intersect_end_pos = vec2row(intersect_end_pos);
    intersect_inds1 = vec2row(intersect_inds1);
    intersect_inds2 = vec2row(intersect_inds2);
end
if(iscolvector(start_pos1) && iscolvector(start_pos2))
    intersect_start_pos = vec2column(intersect_start_pos);
    intersect_end_pos = vec2column(intersect_end_pos);
    intersect_inds1 = vec2column(intersect_inds1);
    intersect_inds2 = vec2column(intersect_inds2);
end

