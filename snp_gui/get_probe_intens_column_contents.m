%function [quartet_num_vec, match_cell, strand_cell] = get_probe_intens_column_contents()
function [quartet_num_vec, match_cell, strand_cell] = get_probe_intens_column_contents()

quartet_num_vec = [1:5];
match_cell = {'PA';'MA';'PB';'MB'};
strand_cell = {'Sense';'Antisense'};