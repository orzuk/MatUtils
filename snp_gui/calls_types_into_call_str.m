%function [call_type_str] = calls_types_into_call_str(call_type_vec)
function [call_type_str] = calls_types_into_call_str(call_type_vec)

[Retention, Non_Inf, No_Call_both, No_Call_first, No_Call_second, LOH, AA_BB_diff, call_type_str_cell] = ...
    call_type_ind_into_call_type_str_ind();

call_type_str = cell(size(call_type_vec));
num_call_types = length(call_type_str_cell);
for i = 1:num_call_types
    call_type = i;
    call_type_ind = find(call_type_vec==call_type);
    if(length(call_type_ind)>0)
        call_type_str = assign_str_to_indices_in_cell(call_type_str, char(call_type_str_cell{i,1}), call_type_ind);
    end
end


