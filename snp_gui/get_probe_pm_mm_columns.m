%function [pm_column] = snp_get_probe_pm_mm_columns(cur_dir, array_name, sample_name)
function [pm_column, mat, pa_vec, ma_vec, pb_vec, mb_vec] = get_probe_pm_mm_columns(cur_dir, array_name, sample_name)

eval(['load(''' cur_dir array_name ''');']);

[pa_vec, ma_vec, pb_vec, mb_vec] = get_probe_pm_mm_ind(column_content_cell, quartet_num_vec, match_cell, strand_cell);
eval(['mat = ' sample_name '_probe_intens_mat;']);
pm_column = reshape(mat(:,[pa_vec pb_vec]), 1, 20*length(mat));