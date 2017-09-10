% transfer index to string 
%function [Retention, Non_Inf, No_Call_both, No_Call_first, No_Call_second, LOH, AA_BB_diff, call_type_str_cell, call_types_vec] = call_type_ind_into_call_type_str_ind()
function [Retention, Non_Inf, No_Call_both, No_Call_first, No_Call_second, ...
    LOH, AA_BB_diff, call_type_str_cell, call_types_vec] = call_type_ind_into_call_type_str_ind()

Retention = 1;
Non_Inf = 2;
No_Call_both = 3;
No_Call_first = 4;
No_Call_second = 5;
LOH = 6;
AA_BB_diff = 7;

call_type_str_cell = {'Retention';'Non-Inf';'NoCall both';'NoCall first'; 'NoCall second';'LOH';'AA/BB-diff'};
call_types_vec = [1:length(call_type_str_cell)];

