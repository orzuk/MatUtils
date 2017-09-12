%function geno_mat = hmm_geno_into_old(in_geno_mat)
function geno_mat = hmm_geno_into_old(in_geno_mat)

[hmm_AA, hmm_AB, hmm_BA, hmm_BB] = genotype_call_into_num_hmm();

[AA, AB, BB, NoCall, call_cell] = genotype_call_into_num();

geno_mat = in_geno_mat;

geno_mat(in_geno_mat==hmm_AA) = AA;
geno_mat(in_geno_mat==hmm_AB) = AB;
geno_mat(in_geno_mat==hmm_BA) = AB;
geno_mat(in_geno_mat==hmm_BB) = BB;