%function [gene_symbols, chr, loc_start, loc_end, exon_start_cell,
%exon_end_cell, gene_descr] = load_gene_loc_exons_ucsc(genome_build)
function [gene_symbols, chr, loc_start, loc_end, exon_start_cell, exon_end_cell, gene_descr] = load_gene_loc_exons_ucsc(genome_build)

load(fullfile('..', 'new_DB', ['refgenes_' genome_build '_with_exons.mat']), 'gene_symbols', 'chr', 'loc_start', 'loc_end', ...
    'exon_start_cell', 'exon_end_cell', 'gene_descr');
