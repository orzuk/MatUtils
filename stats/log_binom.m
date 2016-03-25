% Compute log of binomial coefficient (N over k).
% N and k can be vectors. We use here pre-computed factorial_vec.
% 
% Input: 
% N - number to choose from 
% k - number to choose 
% 
% Output: 
% res - log ( N-choose-k)
% 
function res = log_binom(N, k)

global cumsum_log_vec;

if(length(cumsum_log_vec) < max(N)) % no pre-computation
    res = log_factorial_vec(N) - log_factorial_vec(k) - log_factorial_vec(N-k);
else
    res = cumsum_log_vec(N+1) - cumsum_log_vec(k+1) - cumsum_log_vec(N-k+1); % add one to include log(0) '=' 0
end
