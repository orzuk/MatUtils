% The function loads all the bands which belong to a certain chromsome
function [bands_start_end, bands_length, bands_names] = load_chr_bands(chr)

genome_assembly = get_genome_assembly();

load(['../Database/chr_data_' genome_assembly '.mat'], 'chr_band_names', 'chr_band_start', 'chr_band_end', 'chr_num');

chr_ind = find(chr_num==chr);
num_bands = length(chr_ind);
bands_names = chr_band_names(chr_ind);
bands_start_end = zeros(num_bands, 2);
bands_start_end(:,1) = chr_band_start(chr_ind);
bands_start_end(:,2) = chr_band_end(chr_ind);
bands_length = bands_start_end(:,2)-bands_start_end(:,1);