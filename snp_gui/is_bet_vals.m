%function ret_vec = is_bet_vals(val_vec1, val_vec2, in_val)
function ret_vec = is_bet_vals(val_vec1, val_vec2, in_val)

ret_vec = zeros(length(val_vec1), 1);

bet_ind = find(val_vec1 <= in_val & val_vec2 >= in_val);

ret_vec(bet_ind) = 1;

bet_ind = find(val_vec1 >= in_val & val_vec2 <= in_val);

ret_vec(bet_ind) = 1;