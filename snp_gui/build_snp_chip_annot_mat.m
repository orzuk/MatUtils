%function build_snp_chip_annot_mat(genome_build, chip)
function build_snp_chip_annot_mat(genome_build, chip)

%path(path,'..\src\');

if(strcmp(chip, 'xba')) chip = 'Xba'; end
if(strcmp(chip, 'hind')) chip = 'Hind'; end

f_path = fullfile('..','raw_database');
res_path = fullfile('..','new_DB');


load(fullfile(f_path, [chip '_allele_info_table_' genome_build '.mat']), 'snp_table');
res_mat_f = fullfile(res_path, [chip '_annot_data_' genome_build '.mat']);
first_line = snp_table(1,:);
snp_id_column = strmatch('Probe Set ID', first_line);
rs_column = strmatch('dbSNP RS ID', first_line);
chr_column = strmatch('Chromosome', first_line);
loc_column = strmatch('Physical Position', first_line);
allele_a_column = strmatch('Allele A', first_line);
allele_b_column = strmatch('Allele B', first_line);
strand_column = strmatch(lower('Strand'), lower(first_line));

snp_ids = snp_table(2:end, snp_id_column);
allele_a = snp_table(2:end, allele_a_column);
allele_b = snp_table(2:end, allele_b_column);
chr_vec = snp_table(2:end, chr_column); % change empty indices
chr_loc_vec = snp_table(2:end, loc_column); % change empty indices
rs_ids = snp_table(2:end, rs_column);
strand = snp_table(2:end, strand_column);

chr_vec = assign_empty_str_to_empty_cells(chr_vec);

chr_vec(strcmpi(chr_vec,'x'))=num2cell(23);
chr_vec(strcmpi(chr_vec,'y'))=num2cell(24);
chr_vec(strcmp(chr_vec,'---'))=num2cell(-1);
chr_vec=cell2mat(chr_vec);

%chr_vec = chr_str_cell_into_chr_num_vec(chr_vec);
chr_loc_vec = assign_empty_str_to_empty_cells(chr_loc_vec);

chr_loc_vec(strcmp(chr_loc_vec,'---'))=num2cell(-1);
chr_loc_vec=cell2mat(chr_loc_vec);
%chr_loc_vec = chr_loc_str_cell_into_chr_loc_vec(chr_loc_vec);

[gene_symbols, chr, loc_start, loc_end, exon_start_cell, exon_end_cell, gene_descr] = load_gene_loc_exons_ucsc(genome_build);

% assign for each snp the closest gene
num_snps = length(snp_ids);
snps_gene_symbols = cell(num_snps, 1);
snps_gene_symbols(:) = {''};
%snps_gene_descr = snps_gene_symbols;
snps_gene_dist = snps_gene_symbols;

% for each chromosome find genes close to snps on it
chr_vec_all = [1:24];
num_chr = length(chr_vec_all);
for i = 1:num_chr
    chr_num = chr_vec_all(i);
    genes_chr_ind = find(chr==chr_num);
    num_genes_chr = length(genes_chr_ind);
    snp_chr_ind = find(chr_vec==chr_num);
    num_snp_chr = length(snp_chr_ind);
    dist_mat = zeros(num_snp_chr, num_genes_chr);
    %    dist_cell = cell(num_snp_chr, num_genes_chr);
    if(num_genes_chr & num_snp_chr)
        for j = 1:num_genes_chr
            dist_mat(:,j) = get_dist_from_gene(loc_start(genes_chr_ind(j)), loc_end(genes_chr_ind(j)), ...
                exon_start_cell{genes_chr_ind(j)}, exon_end_cell{genes_chr_ind(j)}, chr_loc_vec(snp_chr_ind));
        end
    end
    % find for each snp its minimal ind
    [temp, min_ind_vec] = min(dist_mat');
    dist_vec = get_mat_spec_ind(dist_mat,[1:size(dist_mat,1)], min_ind_vec);
    dist_cell = gene_dist_vec_into_dist_cell(dist_vec);
    snps_gene_symbols(snp_chr_ind) = gene_symbols(genes_chr_ind(min_ind_vec));
    snps_gene_dist(snp_chr_ind) = dist_cell';
%    snps_gene_descr(snp_chr_ind) = snps_gene_descr_cell;
end

%gene_symbols = snps_gene_symbols;
%gene_descr = snps_gene_descr;

% order according to chr location
ordered_ind = order_by_chr_loc(chr_vec, chr_loc_vec);
snp_ids = snp_ids(ordered_ind);
rs_ids = rs_ids(ordered_ind);
chr_vec = chr_vec(ordered_ind);
chr_loc_vec = chr_loc_vec(ordered_ind);
snps_gene_dist = snps_gene_dist(ordered_ind);
snps_gene_symbols = snps_gene_symbols(ordered_ind);
allele_a = allele_a(ordered_ind);
allele_b = allele_b(ordered_ind);
strand = strand(ordered_ind);

save(res_mat_f, 'snp_ids', 'rs_ids', 'chr_vec', 'chr_loc_vec', 'snps_gene_dist', 'snps_gene_symbols', ...
    'allele_a', 'allele_b', 'genome_build', 'chip', 'strand');
