% Compute log (\sum_i exp(v(i)). 
% Avoid underflows 
% 
function s = my_log_sum_exp(v)

max_v = max(v); 
v = v-max_v; 
s = log(sum(exp(v))) + max_v; 


