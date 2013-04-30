%function [call_type_vec, call_type_str_cell] = calls_ind_into_call_types(no_call_ind1, AB_ind1, AA_ind1, BB_ind1, no_call_ind2, AB_ind2, AA_ind2, BB_ind2, num_snps)
function [call_type_vec, call_type_str_cell] = calls_ind_into_call_types(no_call_ind1, AB_ind1, AA_ind1, BB_ind1, no_call_ind2, AB_ind2, AA_ind2, BB_ind2, num_snps)

[Retention, Non_Inf, No_Call_both, No_Call_first, No_Call_second, LOH, AA_BB_diff, call_type_str_cell] = ...
    call_type_ind_into_call_type_str_ind();
AB_or_no_call_AB = Retention;
same_AA_BB = Non_Inf;
NO_CALL_NO_CALL = No_Call_both;
NO_CALL_else = No_Call_first;
else_NO_CALL = No_Call_second;
diff_AB_AA_or_BB = LOH;
AA_or_BB_diff = AA_BB_diff;


call_type_vec = zeros(1, num_snps);

AA_AA_ind = intersect(AA_ind1, AA_ind2);
BB_BB_ind = intersect(BB_ind1, BB_ind2);
same_AA_BB_ind = [AA_AA_ind; BB_BB_ind];

AA_BB_ind = intersect(AA_ind1, BB_ind2);
AA_AB_ind = intersect(AA_ind1, AB_ind2);
BB_AA_ind = intersect(BB_ind1, AA_ind2);
BB_AB_ind = intersect(BB_ind1, AB_ind2);
AA_or_BB_diff_ind = [AA_BB_ind;AA_AB_ind;BB_AA_ind;BB_AB_ind];

NO_CALL_NO_CALL_ind = intersect(no_call_ind2, no_call_ind1);
NO_CALL_else_ind = no_call_ind1;
[temp, IA, IB] = intersect(NO_CALL_else_ind, no_call_ind2);
NO_CALL_else_ind(IA) = [];
else_NO_CALL_ind = no_call_ind2;
[temp, IA, IB] = intersect(else_NO_CALL_ind, no_call_ind1);
else_NO_CALL_ind(IA) = [];


AB_AB_ind = intersect(AB_ind1, AB_ind2);
no_call_AB_ind = intersect(no_call_ind1, AB_ind2);
AB_or_no_call_AB_ind = [AB_AB_ind; no_call_AB_ind];

AB_AA_ind = intersect(AB_ind1, AA_ind2);
AB_BB_ind = intersect(AB_ind1, BB_ind2);
diff_AB_AA_or_BB_ind = [AB_AA_ind; AB_BB_ind];

call_type_vec(same_AA_BB_ind) = same_AA_BB;
call_type_vec(NO_CALL_else_ind) = NO_CALL_else;
call_type_vec(else_NO_CALL_ind) = else_NO_CALL;
call_type_vec(NO_CALL_NO_CALL_ind) = NO_CALL_NO_CALL;
call_type_vec(AB_or_no_call_AB_ind) = AB_or_no_call_AB; % see this and not no call
call_type_vec(diff_AB_AA_or_BB_ind) = diff_AB_AA_or_BB;
call_type_vec(AA_or_BB_diff_ind) = AA_or_BB_diff;


