%function create_chr_bands_data(genome_build)
function create_chr_bands_data(genome_build, save_path)

chr_num=[];
chr_band_start=[];
chr_band_end=[];
chr_band_names={};

%path = 'E:\Libi\databases\UCSC\';

table = loadCellFile(fullfile('..','raw_database',['cytoBand_' genome_build '.txt']));
chr_num = table(:,1);
% change X and Y chromosomes into 23, 24
x_ind = strmatch('chrX', chr_num);
chr_num(x_ind) = {'chr23'};
x_ind = strmatch('chrx', chr_num);
chr_num(x_ind) = {'chr23'};
y_ind = strmatch('chrY', chr_num);
chr_num(y_ind) = {'chr24'};
y_ind = strmatch('chry', chr_num);
chr_num(y_ind) = {'chr24'};

chr_num = char(chr_num);
chr_num = chr_num(:,4:end);

chr_num = str2num(chr_num);
chr_band_start = cell2mat(table(:,2));
chr_band_end = cell2mat(table(:,3));
chr_band_names = table(:,4);

for i=1:24
    S=char(chr_band_names);
    chr_idx = find(chr_num == i);
    ix=find(S(chr_idx,1)=='q');
    end_p_location(i) = chr_band_end(chr_idx(min(ix)-1));
end

end_p_location = end_p_location';
save(fullfile(save_path, ['chr_data_' genome_build '.mat']), 'chr_num', 'chr_band_start', 'chr_band_end', 'chr_band_names', 'end_p_location');
