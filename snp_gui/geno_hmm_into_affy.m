%function genotype_mat_affy = geno_hmm_into_affy(genotype_mat, data_snp_ids, chip_type, chip_snp_ids_ordered)
function genotype_mat_affy = geno_hmm_into_affy(genotype_mat, data_snp_ids, chip_type, chip_snp_ids_ordered)

genome_assembly = get_genome_assembly();
if strcmp(chip_type,'Hind_and_Xba')
    database_struct = load(['..\database\Hind_annot_data_' genome_assembly '.mat'], 'snp_ids', 'strand');
    database_struct2 = load(['..\database\Xba_annot_data_' genome_assembly '.mat'], 'snp_ids', 'strand');
    
    database_struct.snp_ids=[database_struct.snp_ids; database_struct2.snp_ids];
    database_struct.strand=[database_struct.strand; database_struct2.strand];
    
    [dum idx]= ismember(chip_snp_ids_ordered, database_struct.snp_ids);
    idx=idx(dum);
    
    database_struct.strand=database_struct.strand(idx);
    database_struct.snp_ids=database_struct.snp_ids(idx);
    
else
    database_struct = load(['..\database\' chip_type '_annot_data_' genome_assembly '.mat'], 'snp_ids', 'strand');
end


[C, IA, IB] = intersect_order_by_first_gr(data_snp_ids, database_struct.snp_ids);

strand_sign = zeros(length(C), 1);
strand_sign(strmatch('-', database_struct.strand(IB))) = 1;

genotype_mat_affy = hmm_geno_into_old(genotype_mat);
genotype_mat_affy = geno_into_affy_strand(genotype_mat_affy, strand_sign);
