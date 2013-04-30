%function median_intensity = probe_intens_txt_into_mat(probe_intens_path, array_name, mat_f_name_path, mat_f_name, sample_name)
function median_intensity = probe_intens_txt_into_mat(probe_intens_path, array_name, mat_f_name_path, mat_f_name, sample_name)

save_str = '';
num_header_columns = 2;
num_header_raws = 2;
[quartet_num_vec, match_cell, strand_cell] = get_probe_intens_column_contents();

probe_intens_table = loadCellFile([probe_intens_path array_name]); %%%  '_probe_intens.txt']);
second_line = probe_intens_table(2, :);
snp_id_column = strmatch('Probe Set', second_line);
sample_snp_id = probe_intens_table(3:end, snp_id_column);
snp_id = sample_snp_id;
num_snp = length(snp_id);

probe_intens_table_mat = probe_intens_table(num_header_raws+1:end,num_header_columns+1:end);
null_ind = strcmp('null', probe_intens_table_mat);
probe_intens_table_mat(null_ind) = {'-1'};
probe_intens_table_mat = cell_of_nums_to_mat(probe_intens_table_mat);
% extract the intensities
% 5 quartet for each strand
num_columns = length(quartet_num_vec)*length(match_cell)*length(strand_cell);

column_ind = 1;
sample_probe_intens_mat = zeros(num_snp, num_columns);
for j = 1:length(quartet_num_vec)
    for t = 1:length(match_cell)
        for q = 1:length(strand_cell)
            column_name = [sample_name '_Quartet' num2str(quartet_num_vec(j)) ' - ' ...
                char(match_cell{t}) '(' char(strand_cell{q}) ')' ]; %95Ax_Quartet3 - PB(Sense)
            column_loc = strmatch(column_name, second_line);
            sample_probe_intens_mat(:, column_ind) = ...
                probe_intens_table_mat(:, column_loc - num_header_columns);
            column_content_cell{column_ind} = [j t q];
            column_ind = column_ind+1;
        end
    end
end

[pa_vec, ma_vec, pb_vec, mb_vec] = get_probe_pm_mm_ind(column_content_cell, quartet_num_vec, match_cell, strand_cell);
pm_column = reshape(sample_probe_intens_mat(:,[pa_vec pb_vec]), 1, length([pa_vec pb_vec])*length(sample_probe_intens_mat));
median_intensity = median(pm_column);

eval([ sample_name '_probe_intens_mat = sample_probe_intens_mat;']);
save_str = [save_str ',''' sample_name '_probe_intens_mat'''];
save_str = [save_str ',''' 'column_content_cell'''];

save_str = [save_str ',''' 'quartet_num_vec'''];
save_str = [save_str ',''' 'match_cell'''];
save_str = [save_str ',''' 'strand_cell'''];

% save into matlab file
eval(['save(''' mat_f_name_path  '\' mat_f_name '''' save_str ');']);



