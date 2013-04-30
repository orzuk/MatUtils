% Compute log-factorial for a vector. 
% Not only one number like in factorial of
% matlab. Can we use here steirling ?
%
% Input: 
% N - vector of integers
% 
% Output: 
% real_res_vec - log factorial 
%
function[real_res_vec]= log_factorial_vec(N)

global cumsum_log_vec; % why global? 

[unisorted_N ind rev_ind] = unique(N); % Take all input numbers

arr = zeros(1,length(unisorted_N)); % Array with the results
% arr(1) = log(factorial(unisorted_N(1)))  - overflow problems !!!
arr(1) = cumsum_log_vec(unisorted_N(1)); % arr(1) = sum(log(1:unisorted_N(1)));
for i=2:length(unisorted_N) % loop over all numbers
    arr(i) = cumsum_log_vec(unisorted_N(i)); 
    % arr(i) = arr(i-1) + sum(log(unisorted_N(i-1)+1:unisorted_N(i)));
end
res_vec = arr(rev_ind); % fill back to vector 
real_res_vec = res_vec; % Change dimensions if all are the same
if(size(real_res_vec, 1) > 1) % transpose 
    real_res_vec = real_res_vec';
end
