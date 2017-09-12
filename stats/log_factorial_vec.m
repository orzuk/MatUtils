% Compute log-factorial for a vector. 
% Not only one number like in factorial of Matlab. Can we use here steirling ?
%
% Input: 
% N - vector of integers
% 
% Output: 
% real_res_vec - log (N!) 
%
function real_res_vec = log_factorial_vec(N)

global cumsum_log_vec; % why global? 

if(length(cumsum_log_vec) == 1) % cumsum_log_vec not initialized yet
    cumsum_log_vec = cumsum([0 log(1:max(N))]);
end
if(length(cumsum_log_vec) < max(N)) % cumsum_log_vec initialized but too short
    cumsum_log_vec = cumsum([0 log(1:max(N))]); % can make more efficient: compute only new values !! 
end

[unisorted_N, ~, rev_ind] = unique(N); % Take all input numbers

arr = zeros(1,length(unisorted_N)); % Array with the results
% arr(1) = log(factorial(unisorted_N(1)))  - overflow problems !!!
%arr(1) = cumsum_log_vec(unisorted_N(1)); % arr(1) = sum(log(1:unisorted_N(1)));
for i=1:length(unisorted_N) % loop over all numbers
    arr(i) = cumsum_log_vec(unisorted_N(i)+1); 
    % arr(i) = arr(i-1) + sum(log(unisorted_N(i-1)+1:unisorted_N(i)));
end
real_res_vec = arr(rev_ind); % fill back to vector.  Change dimensions if all are the same
if(size(real_res_vec, 1) > 1) % transpose 
    real_res_vec = real_res_vec';
end
