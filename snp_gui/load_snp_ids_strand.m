%function [ret_strand ret_snp_ids] = load_snp_ids_strand(in_snp_ids, chip)
function [ret_strand ret_snp_ids] = load_snp_ids_strand(in_snp_ids, chip)


AssignAllGlobalConstants;
num_snps = length(in_snp_ids);

ret_strand = cell(num_snps, 1);
ret_strand(:) = {''};

load(['..\database\' chip '_annot_data_' genome_assembly '.mat'], 'snp_ids', 'strand');

ret_snp_ids = snp_ids;
if(num_snps>0)
    [C, IA, IB] = intersect(in_snp_ids, snp_ids);

    ret_strand(IA) = strand(IB);
    ret_snp_ids = in_snp_ids;
end