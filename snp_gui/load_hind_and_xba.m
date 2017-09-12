function chip_annot=load_Hind_and_Xba

genome_assembly = get_genome_assembly();

Hind_chip_annot=load(fullfile('..','database',['Hind_annot_data_' genome_assembly '.mat']));
Xba_chip_annot=load(fullfile('..','database',['Xba_annot_data_' genome_assembly '.mat']));

chip_annot.genome_build = Hind_chip_annot.genome_build;
chip_annot.chip= 'Hind_and_Xba';
chip_annot.snp_ids = [Hind_chip_annot.snp_ids; Xba_chip_annot.snp_ids];
chip_annot.allele_a = [Hind_chip_annot.allele_a; Xba_chip_annot.allele_a];
chip_annot.allele_b = [Hind_chip_annot.allele_b; Xba_chip_annot.allele_b];
chip_annot.rs_ids = [Hind_chip_annot.rs_ids; Xba_chip_annot.rs_ids];
chip_annot.strand = [Hind_chip_annot.strand; Xba_chip_annot.strand];
chip_annot.chr_vec = [Hind_chip_annot.chr_vec; Xba_chip_annot.chr_vec];
chip_annot.chr_loc_vec = [Hind_chip_annot.chr_loc_vec; Xba_chip_annot.chr_loc_vec];
chip_annot.snps_gene_symbols = [Hind_chip_annot.snps_gene_symbols; Xba_chip_annot.snps_gene_symbols];
chip_annot.snps_gene_dist = [Hind_chip_annot.snps_gene_dist; Xba_chip_annot.snps_gene_dist];

chr_loc_idx=[chip_annot.chr_loc_vec chip_annot.chr_vec [1:length(chip_annot.chr_loc_vec)]'];
chr_loc_idx=sortrows(chr_loc_idx,[2 1]); %sort by chr and then by loc
idx=chr_loc_idx(:,3);

chip_annot.snp_ids = chip_annot.snp_ids(idx);
chip_annot.allele_a = chip_annot.allele_a(idx);
chip_annot.allele_b = chip_annot.allele_b(idx);
chip_annot.rs_ids = chip_annot.rs_ids(idx);
chip_annot.strand = chip_annot.strand(idx);
chip_annot.chr_vec = chip_annot.chr_vec(idx);
chip_annot.chr_loc_vec = chip_annot.chr_loc_vec(idx);
chip_annot.snps_gene_symbols = chip_annot.snps_gene_symbols(idx);
chip_annot.snps_gene_dist = chip_annot.snps_gene_dist(idx);







