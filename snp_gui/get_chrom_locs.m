% This function outputs the chromosomes and chromosomal locations of input snps. 
% It also outputs the pres_ind vec, which says what SNPs are present on the
% chip (some snps might be missing)
%function [chr, locs, pres_ind] = get_chrom_locs(snp_ids, chip)
function [chr, locs, pres_ind] = get_chrom_locs(snp_ids, chip)

genome_assembly = get_genome_assembly();
chip_annot=load(fullfile('..','database',[chip '_annot_data_' genome_assembly '.mat']), 'snp_ids', 'chr_vec', 'chr_loc_vec');

[C, pres_ind, IB] = intersect_order_by_first_gr(snp_ids, chip_annot.snp_ids);

chr = chip_annot.chr_vec(IB);
locs = chip_annot.chr_loc_vec(IB);