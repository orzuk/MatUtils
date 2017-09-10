%function ordered_ind = order_by_chr_loc(chr_num_vec, chr_loc_vec)
function ordered_ind = order_by_chr_loc(chr_num_vec, chr_loc_vec)
% order data by chr and location
[temp, chr_ind_sorted] = sort(chr_num_vec);
ordered_ind = chr_ind_sorted;
chr_vec = unique(chr_num_vec);
for i = 1:length(chr_vec)
    chr_ind = find(chr_num_vec(chr_ind_sorted)==chr_vec(i));
    [temp, loc_in_chr_sorted_ind] = sort(chr_loc_vec(chr_ind_sorted(chr_ind)));
    ordered_ind(chr_ind) = ordered_ind(chr_ind(loc_in_chr_sorted_ind)) ;
end
