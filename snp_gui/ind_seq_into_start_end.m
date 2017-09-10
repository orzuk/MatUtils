%function [start_vec, end_vec] = ind_seq_into_start_end(ind_seq)
function [start_vec, end_vec] = ind_seq_into_start_end(ind_seq)

% find where ind is not consecutive

diff_vec = diff(ind_seq);
diff_ind = find(diff_vec ~= 1);
num_diff = length(diff_ind);
start_vec = zeros(num_diff+1,1);
end_vec = zeros(num_diff+1,1);
if(num_diff>0)
    curr_start_ind = 1;
    curr_end_ind = 1;
    for i = 1:num_diff
        curr_end_ind = diff_ind(i);
        start_vec(i) = ind_seq(curr_start_ind);
        end_vec(i) = ind_seq(curr_end_ind);
        curr_start_ind = curr_end_ind+1;
    end
    start_vec(end) = ind_seq(curr_start_ind);
    end_vec(end) = ind_seq(end);
else
    start_vec = [ind_seq(1)];
    end_vec = [ind_seq(end)];
end