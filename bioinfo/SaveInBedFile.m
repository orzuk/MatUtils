% Save data in .bed format:
% chr, start, end, gene-name, score, strand, thickstart, thickend, RGB, #exons, exon-lengths (relative), exon-start,
%
function SaveInBedFile(output_bed_file_name, chr_vec, pos_start_vec, pos_end_vec, ...
    genome_version, gene_names, strand, ...
    scores, num_exons, exon_starts, exon_lengths)


num_lines = length(chr_vec);
%data_struct = cell(num_lines, 6); % gene-name, strand, score, #exons, exon-start, exon-lengths
%data_struct(:,2) = num2strand(strand); % write '+' or '-' strand

data_struct = [gene_names mat2cell(scores, ones(num_lines,1)) ...
    mat2cell(num2strand(strand), ones(num_lines,1)), ...
    mat2cell(pos_start_vec, ones(num_lines,1)), ...
    mat2cell(pos_end_vec, ones(num_lines,1)), ...
    mat2cell(repmat('0,0,0', num_lines, 1), ones(num_lines,1)), ...
    mat2cell(num_exons, ones(num_lines,1))];
for i=1:num_lines
%     do_line = i
%     exon_lengths_are = exon_lengths{i}'
    data_struct{i,8} = cell2vec(num2str_cell(exon_lengths{i}'), ',');
    data_struct{i,9} = cell2vec(num2str_cell(exon_starts{i}'), ',');
end

add_chr_flag = 1;
save_regions_mat_as_text(output_bed_file_name, chr_vec, pos_start_vec, pos_end_vec, ...
    genome_version, data_struct, add_chr_flag);

