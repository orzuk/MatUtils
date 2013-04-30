%function chr_num_vec = chr_str_cell_into_chr_num_vec(chr_str_cell)
function chr_num_vec = chr_str_cell_into_chr_num_vec(chr_str_cell)

chr_num_vec = chr_str_cell;
chr_num_vec  = change_strs_in_cell(chr_num_vec , {'X';'x'}, '23');
chr_num_vec  = change_strs_in_cell(chr_num_vec , {'Y';'y'}, '24');
chr_num_vec  = change_strs_in_cell(chr_num_vec , {''}, '-1');
chr_num_vec = cell2mat(chr_num_vec);
