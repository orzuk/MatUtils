%function [call_ind_str] = calls_ind_into_call_str(call_ind_vec)
function [call_ind_str] = calls_ind_into_call_str(call_ind_vec)

[AA, AB, BB, NoCall, call_cell] = genotype_call_into_num();

call_ind_str = cell(size(call_ind_vec));
num_call_ind = length(call_cell);
for i = 1:num_call_ind
    call = i;
    call_ind = find(call_ind_vec==call);
    if(length(call_ind)>0)
        call_ind_str = assign_str_to_indices_in_cell(call_ind_str, char(call_cell{i,1}), call_ind);
    end
end


