% Here we check that in the BSC everything is ok : 

xi = [1,1]';
U = [0.5 0.5; 0.5 0.5];
T = [1 -1; -1 1];
delta = 0.2;

M = U + delta * T;

eps = 0.1; R = [1-eps eps; eps 1-eps];

pi = [0.5 0.5];

BSC_KL_first_order_A_M = (1/2) .* [xi' * ((R' * T * R) .* (log((1/2)*R'*U*R) ))] * xi
BSC_KL_first_order_HIGH_SNR =  xi' * [ diag(log(pi)) * T' * diag(pi) * M - (diag(pi)*M*T+T'*diag(pi)*M) .* log(diag(pi)*M)       ] * xi


p = 0.5-delta;
BSC_ORDER = 2 * (1-2*p) * log((1-p)/p)