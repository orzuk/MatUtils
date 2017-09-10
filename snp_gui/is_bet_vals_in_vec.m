%function ret_vec = is_bet_vals_in_vec(val1, val2, in_val_vec)
function ret_vec = is_bet_vals_in_vec(val1, val2, in_val_vec)

ret_vec = zeros(length(in_val_vec), 1);

if(val1>val2)
    temp = val1;
    val1 = val2;
    val2 = temp;
end
bet_ind = find(in_val_vec >= val1 & in_val_vec <= val2);

ret_vec(bet_ind) = 1;
