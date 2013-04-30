% Compute log of binomial coefficient N over M.
% N and M can be vectors. We use here factorial_vec.
function[res]= log_binom(N, M)

global cumsum_log_vec;

%res = log_factorial_vec(N) - log_factorial_vec(M) - log_factorial_vec(N-M);
res = cumsum_log_vec(N+1) - cumsum_log_vec(M+1) - cumsum_log_vec(N-M+1); % add one to include log(0) '=' 0
