% Just like ReadWigFileToMat, with another name
function [chr_vec pos_start_vec pos_end_vec data] = ...
    ReadBedFileToMat(wig_file, genome_version, mat_outfile, numeric_flag, varargin)

if(~exist('mat_outfile', 'var'))
    mat_outfile = [];
end
if(~exist('numeric_flag', 'var'))
    numeric_flag = [];
end

old_wig = 0;
[chr_vec pos_start_vec pos_end_vec data] = ...
    ReadWigFileToMat(wig_file, genome_version, mat_outfile, numeric_flag, varargin);
if(~old_wig) % New: here read a full bed format file
    data2 = data;
    data.gene_names = data2(:,1); % fourth field
    data.score = cell2mat(str2num_cell(data2(:,2)));
    data.strand = data2(:,3);
    data.thick_starts = cell2mat(str2num_cell(data2(:,4)));
    data.thick_ends = cell2mat(str2num_cell(data2(:,5)));
    data.RGB = str2nums_cell(data2(:,6));
    data.num_exons = cell2mat(str2num_cell(data2(:,7)));
    data.exon_lengths = str2nums_cell(data2(:,8));
    data.exon_starts = str2nums_cell(data2(:,9));
    if(~isempty(mat_outfile))
        save(mat_outfile, 'data', '-append');
    end
end

