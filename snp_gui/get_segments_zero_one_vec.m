%function [joined_ones_vec] = get_segments_zero_one_vec(zero_one_vec)
function [joined_ones_vec] = get_segments_zero_one_vec(zero_one_vec)

if(size(zero_one_vec,1) ~=1) zero_one_vec = zero_one_vec'; end
vec_len = length(zero_one_vec);
zero_one_vec = [0 zero_one_vec 0];
Segs = zero_one_vec(1:end-1) - zero_one_vec(2:end);
SegsStarts = find(Segs == -1); SegsEnds = find(Segs == 1)-1;
num_segs = length(SegsStarts);
joined_ones_vec = zeros(2,num_segs);
joined_ones_vec(1,:) = SegsEnds-SegsStarts+1;
joined_ones_vec(2,:) = SegsStarts;
stretch_len_vec = zeros(1, vec_len);

