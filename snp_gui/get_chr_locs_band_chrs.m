%function bands_name_cell = get_chr_locs_band_chrs(chr_vec, chr_loc_vec)
function bands_name_cell = get_chr_locs_band_chrs(chr_vec, chr_loc_vec)

bands_name_cell = cell(length(chr_vec),1);
bands_name_cell(:) = {''};
num_chr = 24;
for i = 1:num_chr
    chr_ind = find(chr_vec==i);
    if(length(chr_ind)>0)
        [bands_start_end, bands_length, bands_names] = load_chr_bands(i);
        chr_bands_name_cell = get_chr_locs_band(chr_loc_vec(chr_ind), bands_start_end, bands_names);
        chr_bands_name_cell = add_pref_to_cell_str(chr_bands_name_cell, num2str(i));
        bands_name_cell(chr_ind) = chr_bands_name_cell;
    end
end
