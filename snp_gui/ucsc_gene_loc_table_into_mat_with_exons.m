%function ucsc_gene_loc_table_into_mat_with_exons_script(genome_build)
function ucsc_gene_loc_table_into_mat_with_exons(genome_build)

%path(path,'..\src\');
load(fullfile('..','raw_database',['refgenes_' genome_build '_cell.mat']), 'gene_table');
%res_path = '..\new_DB\';

nm_ind = strmatch('NM', gene_table(:, 1));
all_gene_symbols = gene_table(nm_ind,11);
all_gene_chr = gene_table(nm_ind,2);
all_gene_start = gene_table(nm_ind,6);
all_gene_end = gene_table(nm_ind,7);
all_gene_num_exon = gene_table(nm_ind,8);
all_gene_start_exon = gene_table(nm_ind,9);
all_gene_end_exon = gene_table(nm_ind,10);
all_gene_descr = gene_table(nm_ind,12);


for i = 1:length(all_gene_symbols) 
    if(isempty(all_gene_symbols{i})) 
        all_gene_symbols{i} = '';
    elseif ~ischar(all_gene_symbols{i})
         all_gene_symbols{i}=num2str(all_gene_symbols{i});
    end
end
[gene_symbols,I,J]  = unique(all_gene_symbols);
gene_chr = all_gene_chr(I);
gene_start = all_gene_start(I);
gene_end = all_gene_end(I);
gene_num_exon = all_gene_num_exon(I);
gene_start_exon = all_gene_start_exon(I);
gene_end_exon = all_gene_end_exon(I);
gene_descr = all_gene_descr(I);

% delete empty gene
t = strcmp('', gene_symbols);
t = find(t==1);
gene_symbols(t) = [];
gene_chr(t) = [];
gene_start(t) = [];
gene_end(t) = [];
gene_num_exon(t) = [];
gene_start_exon(t) = [];
gene_end_exon(t) = [];
gene_descr(t) = [];


num_genes = length(gene_chr);
chr = zeros(num_genes, 1);
loc_start = cell2mat(gene_start);
loc_end = cell2mat(gene_end);
num_exons = cell2mat(gene_num_exon);

for i = 1:num_genes
    gene_chr_str = char(gene_chr{i});
    gene_chr_str = gene_chr_str(4:end);
    t = strfind(gene_chr_str, '_');
    if(length(t)>0)
        gene_chr_str = gene_chr_str(1:min(t)-1);
    end
    gene_chr{i} = gene_chr_str;
end

chr = chr_str_cell_into_chr_num_vec(gene_chr);
exon_start_cell = cell(num_genes, 1);
exon_end_cell = cell(num_genes, 1);
% now extract exons start-end
for i = 1:num_genes
    exon_start_str = char(gene_start_exon{i});
    exon_end_str = char(gene_end_exon{i}); 
    %look like this: 4692,4901,5810,6631,6918,7231,7605,7924,8242,14754,
    b = strfind(exon_start_str,'"');
    exon_start_str(b) = ' ';
    b = strfind(exon_start_str,',');
    exon_start_str(b) = ' ';
    b = strfind(exon_end_str,'"');
    exon_end_str(b) = ' ';
    b = strfind(exon_end_str,',');
    exon_end_str(b) = ' ';
    words_start = parse_words_from_line(exon_start_str);
    exon_start_cell{i} = cell_of_nums_to_mat(words_start);
    words_end = parse_words_from_line(exon_end_str);
    exon_end_cell{i} = cell_of_nums_to_mat(words_end);
end

% order according to chr location
ordered_ind = order_by_chr_loc(chr, loc_start);

gene_symbols = gene_symbols(ordered_ind);
chr = chr(ordered_ind);
loc_start = loc_start(ordered_ind);
loc_end = loc_end(ordered_ind);
gene_descr = gene_descr(ordered_ind);
exon_start_cell = exon_start_cell(ordered_ind);
exon_end_cell = exon_end_cell(ordered_ind);
save(fullfile('..','new_DB', ['refgenes_' genome_build '_with_exons.mat']), 'gene_symbols', 'chr', 'loc_start', 'loc_end', 'gene_descr', ...
    'exon_start_cell', 'exon_end_cell');
