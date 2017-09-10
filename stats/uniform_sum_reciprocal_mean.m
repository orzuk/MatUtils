% Compute the mean of 1/(sum_{i} p_i)
% 
% Input: 
% max_n - number of uniform r.v.s.
% iters - num. of iterations to simulate
% 
% Outuput: 
% mu - mean of  1/(sum_{i} p_i)
% 
function mu = uniform_sum_reciprocal_mean(max_n, iters)

mu = zeros(max_n, 1); 

for i=1:iters
    r = rand(max_n, 1); % sample from uniform distribution 
    mu = mu + 1 ./ cumsum(r);
end
mu = mu ./ iters;


