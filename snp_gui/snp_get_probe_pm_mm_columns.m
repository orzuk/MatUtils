%function [pa_vec, ma_vec, pb_vec, mb_vec] = snp_get_probe_pm_mm_columns(proj_name)
function [pa_vec, ma_vec, pb_vec, mb_vec] = snp_get_probe_pm_mm_columns(proj_name)

[column_content_cell, quartet_num_vec, match_cell, strand_cell] = load_snp_probe_intens_info(proj_name);

pa_ind = strmatch('PA',match_cell);
ma_ind = strmatch('MA',match_cell);
pb_ind = strmatch('PB',match_cell);
mb_ind = strmatch('MB',match_cell);
%match_cell = {'PA';'MA';'PB';'MB'};
num_quartet = length(quartet_num_vec);
num_strands = length(strand_cell);
match_cell = cell(4,num_strands);
match_cell(:) = {zeros(1, num_quartet)};

num_columns = length(column_content_cell);
for i = 1:num_columns
    content_vec = column_content_cell{i};
    temp = match_cell{content_vec(2), content_vec(3)};
    temp(content_vec(1)) = i;
    match_cell{content_vec(2), content_vec(3)} = temp;
end

pa_vec = [];
ma_vec = [];
pb_vec = [];
mb_vec = [];
for i = 1:num_strands
    pa_vec = [pa_vec match_cell{pa_ind,i}];
    ma_vec = [ma_vec match_cell{ma_ind,i}];
    pb_vec = [pb_vec match_cell{pb_ind,i}];
    mb_vec = [mb_vec match_cell{mb_ind,i}];
end
