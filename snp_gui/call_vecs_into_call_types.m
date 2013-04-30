%function [call_type_vec, call_type_str_cell] = call_vecs_into_call_types(calls_vec1, calls_vec2)
function [call_type_vec, call_type_str_cell] = call_vecs_into_call_types(calls_vec1, calls_vec2)



num_snps = length(calls_vec1);

[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(calls_vec1);
[no_call_ind2, AB_ind2, AA_ind2, BB_ind2] = calls_vec_into_ind(calls_vec2);

[call_type_vec, call_type_str_cell] = calls_ind_into_call_types...
    (no_call_ind1, AB_ind1, AA_ind1, BB_ind1, no_call_ind2, AB_ind2, AA_ind2, BB_ind2, num_snps);
call_type_vec = call_type_vec';
