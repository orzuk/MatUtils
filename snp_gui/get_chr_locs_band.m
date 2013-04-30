% This function gets as input a vector of chromosomal locations, and the
% start and end of each band. It outputs the correct band of each location
%function bands_name_cell = get_chr_locs_band(chr_loc_vec, bands_start_end, bands_names)
function bands_name_cell = get_chr_locs_band(chr_loc_vec, bands_start_end, bands_names)

if(size(chr_loc_vec, 1)~=1)
    chr_loc_vec = chr_loc_vec';
end
bands_start = bands_start_end(:,1); 
bands_end = bands_start_end(:,2); 

% find start loc that are before the input chr location
start_ind = find(repmat(bands_start, 1, length(chr_loc_vec)) <= repmat(chr_loc_vec, length(bands_start), 1));
start_mat = zeros(length(bands_start), length(chr_loc_vec));
start_mat(start_ind) = 1;

% find end loc that are after the input chr location
end_ind = find(repmat(bands_end, 1, length(chr_loc_vec)) >= repmat(chr_loc_vec, length(bands_end), 1));
end_mat = zeros(length(bands_end), length(chr_loc_vec));
end_mat(end_ind) = 1;

band_mat = start_mat & end_mat;
[band_ind_vec, j] = find(band_mat==1);

bands_name_cell = cell(length(chr_loc_vec), 1); bands_name_cell(:) = {''};
bands_name_cell(j, 1) = bands_names(band_ind_vec,:);