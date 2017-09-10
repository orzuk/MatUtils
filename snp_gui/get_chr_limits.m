% The function gets the limits of a chromosome, from zero
% to the length of the chromosome
function [xmin, xmax] =  get_chr_limits(chr_num)

[bands_start_end, bands_length, bands_names] = load_chr_bands(chr_num);

xmin = min(min(bands_start_end));
xmax = max(max(bands_start_end));