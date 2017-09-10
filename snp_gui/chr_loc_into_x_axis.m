% The function converts the chromosomal locations into actuall x-axis which
% is used for plotting
function [plot_location, chr_start_loc, end_p_location, ind] = chr_loc_into_x_axis(chr_num_vec, chr_loc_vec, chr_vec)

num_chr = length(chr_vec);
end_p_location = load_end_p_loc();
end_p_location = end_p_location(chr_vec);
% work only with chr. from chr_vec
[B,I,J] = unique(chr_num_vec);
ind = zeros(1, length(chr_num_vec));
[C, IA, IB] = intersect(B, chr_vec);
for i = 1:length(C) % add indices of this chr.
    chr_B_ind = IA(i);
    ind_to_add = find(J==chr_B_ind);
    ind(ind_to_add) = 1;
end
chr_loc_vec(find(ind==0)) = [];
chr_num_vec(find(ind==0)) = [];
[xmin_vec, xmax_vec] =  get_chrs_limits(chr_vec);
plot_location = chr_loc_vec;

chr_change = find(diff(chr_num_vec) ~= 0);
chr_start_loc = ones(1, num_chr);
chr_end_loc = zeros(1, num_chr);
chr_end_loc(1) = xmax_vec(chr_vec(1));
chr_len_vec = xmax_vec(chr_vec);
for i=1:length(chr_change)-1
    chr_end_loc(i+1) = chr_end_loc(i)+ chr_len_vec(i+1);
    end_p_location(i+1) = end_p_location(i+1)+ chr_end_loc(i);
    chr_start_loc(i+1) = chr_start_loc(i+1)+ chr_end_loc(i);
    plot_location(chr_change(i)+1:chr_change(i+1)) = plot_location(chr_change(i)+1:chr_change(i+1)) + chr_end_loc(i);
end

chr_end_loc(end) = chr_end_loc(end-1)+ chr_len_vec(end);
plot_location(chr_change(end)+1:end) = plot_location(chr_change(end)+1:end) + chr_end_loc(end-1);
chr_start_loc(end) = chr_start_loc(end)+ chr_end_loc(end-1);
end_p_location(end) = end_p_location(end)+ chr_end_loc(end-1);



