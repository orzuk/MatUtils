%function ret_geno = geno_into_affy_strand(geno_mat, strand)
function ret_geno = geno_into_affy_strand(geno_mat, strand)

[AA, AB, BB, NoCall, call_cell] = genotype_call_into_num();

if(size(strand, 1)==1) strand = strand'; end
strand_mat = repmat(strand, 1, size(geno_mat, 2));
ret_geno = geno_mat;

ret_geno(find(geno_mat == AA & strand_mat==1)) = BB;
ret_geno(find(geno_mat == BB & strand_mat==1)) = AA;