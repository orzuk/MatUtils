%function [xmin_vec, xmax_vec] =  get_chrs_limits(chr_num_vec)
function [xmin_vec, xmax_vec] =  get_chrs_limits(chr_num_vec)

num_chr = length(chr_num_vec);
for i = 1:num_chr
    [xmin_vec(i), xmax_vec(i)] =  get_chr_limits(chr_num_vec(i));
end
