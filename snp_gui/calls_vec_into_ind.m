%function [no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(calls_vec)
function [no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(calls_vec)

[AA, AB, BB, NoCall] = genotype_call_into_num();

no_call_ind1 = find(calls_vec == NoCall);
AB_ind1 = find(calls_vec == AB);
AA_ind1 = find(calls_vec == AA);
BB_ind1 = find(calls_vec == BB);


