% Compute complement of a set of intervals
function [comp_start_pos comp_end_pos] = ... % take complement of second set
    intervals_complement(start_pos, end_pos, min_val, max_val, is_sorted, varargin)

[start_pos end_pos] = intervals_merge(start_pos, end_pos); 

if(~exist('is_sorted', 'var') || isempty(is_sorted))
    is_sorted = 0; 
end
if(~is_sorted) % always do sorting in matlab (faster)
    tmp_max = max(start_pos, end_pos); tmp_min = min(start_pos, end_pos); start_pos = tmp_min; end_pos = tmp_max;
    [start_pos sort_perm] = sort(start_pos);
    end_pos = end_pos(sort_perm); % sort positions vec by start
end

comp_start_pos = end_pos(1:end-1);
comp_end_pos = start_pos(2:end);
if(start_pos(1) > min_val)
    is_row_flag = isrow(comp_start_pos); % use matlab's isrow (no need for isrowvector) 
    comp_start_pos = [vec2row(min_val) vec2row(comp_start_pos)];
    comp_end_pos = [start_pos(1) vec2row(comp_end_pos)];
    if(~is_row_flag)
        comp_start_pos = comp_start_pos'; 
        comp_end_pos = comp_end_pos'; 
    end
end
if(end_pos(end) < max_val)
    comp_start_pos(end+1) = end_pos(end);
    comp_end_pos(end+1) = max_val;
end

pos_len_inds = find(comp_end_pos > comp_start_pos); % Remove intervals of length 0
comp_start_pos = comp_start_pos(pos_len_inds); 
comp_end_pos = comp_end_pos(pos_len_inds); 


