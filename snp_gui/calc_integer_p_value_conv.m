%function p_vals_up = calc_integer_p_value_conv(int_mat)
function p_vals_up = calc_integer_p_value_conv(int_mat)

time_start = cputime;

%figure; imagesc(int_mat); colorbar; title('not-permuted');


min_val = min(min(int_mat));
max_val = max(max(int_mat));
val_vec = [min_val:max_val];
num_vals = max_val-min_val+1;
num_rows = size(int_mat,1);
num_columns = size(int_mat,2);
vals_p_values = zeros(num_vals, num_columns);
for i = 1:num_columns
    temp = hist(int_mat(:,i), val_vec);
    vals_p_values(:, i) = temp'/num_rows; % for each column keeps the fraction of rows in which each val appears
end

sum_vec = sum(int_mat')'; % sum of each row
p_vals_vec = vals_p_values(:,1);
for i = 2:num_columns
    p_vals_vec = conv(p_vals_vec, vals_p_values(:,i)); % calculates the probability to get each possible value as the sum of a row when permuting each column
end
num_poss_vals = length(p_vals_vec);
poss_vals_min = min_val*num_columns;
poss_vals_max = poss_vals_min+num_poss_vals-1;
poss_vals_vec = [poss_vals_min:poss_vals_max];
poss_vals_vec_p_vals = p_vals_vec;
poss_vals_vec_p_vals_cum = zeros(1, num_poss_vals);

poss_vals_vec_p_vals_cum = cumsum(poss_vals_vec_p_vals);
if(size(poss_vals_vec_p_vals_cum,1) ~=1) poss_vals_vec_p_vals_cum = poss_vals_vec_p_vals_cum'; end
poss_vals_vec_p_vals_cum = [1 1-poss_vals_vec_p_vals_cum(1:end-1)]; % the chance to get 
% poss_vals_vec(i) or more = 1-poss_vals_vec_p_vals_cum(i-1)

% poss_vals_vec_p_vals_cum(1) = 1;
% 
% for i = 2:num_poss_vals
%     poss_vals_vec_p_vals_cum(i) = 1-sum(poss_vals_vec_p_vals(1:i-1));
% end
[unique_sum_vec, I, J] = unique(sum_vec); % sum_vec = unique_sum_vec(J)
p_vals_up = ones(size(unique_sum_vec));
[C, IA, IB] = intersect(unique_sum_vec, poss_vals_vec);
p_vals_up(IA) = poss_vals_vec_p_vals_cum(IB);
p_vals_up = p_vals_up(J);
% 
% for i = 1:num_poss_vals
%     p_vals_up(find(sum_vec==poss_vals_vec(i))) = poss_vals_vec_p_vals_cum(i);
% end

time_end = cputime;
time = time_end-time_start;
