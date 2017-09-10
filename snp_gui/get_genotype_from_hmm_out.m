% The function extracts the (discrete) geneoypes from the HMM output probabilities
function genotype_vec = get_genotype_from_hmm_out(VitStruct, ProbsStruct);

[AA, AB, BB, NoCall, call_cell] = genotype_call_into_num();

num_snp = length(VitStruct.alpha_genotype);

% VitStruct.alpha_genotype: zero for A, one for B
geno_sum_vec = VitStruct.alpha_genotype+VitStruct.beta_genotype;
genotype_vec = zeros(num_snp, 1);
genotype_vec(find(geno_sum_vec==0)) = AA;
genotype_vec(find(geno_sum_vec==1)) = AB;
genotype_vec(find(geno_sum_vec==2)) = BB;