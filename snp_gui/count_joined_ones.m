%function [joined_ones_vec, stretch_len_vec = count_joined_ones(zero_one_vec)
% joined_ones_vec: first line lists the lengths of the stretches, second
% line the location of the beginning of the stretches
function [joined_ones_vec, stretch_len_vec] = count_joined_ones(zero_one_vec)

vec_len = length(zero_one_vec);
joined_ones_vec = get_segments_zero_one_vec(zero_one_vec);
stretch_len_vec = zeros(1, vec_len);
num_segs = size(joined_ones_vec, 2);
for i = 1:num_segs
    start_seg = joined_ones_vec(2,i);
    seg_len = joined_ones_vec(1,i);
    stretch_len_vec(start_seg:start_seg+seg_len-1) = seg_len;
end
