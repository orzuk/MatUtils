% Try Ido's idea with transfer-matrix for HMM 
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary
p = 0.1; eps = 0.014;

A = [   (p^2 + (1-p)^2)*(eps^2 + (1-eps)^2)    2*p*(1-p)*(eps^2 + (1-eps)^2); 2*p*(1-p)*2*eps*(1-eps)  (p^2 + (1-p)^2)*2*eps*(1-eps)]


lambda = eig(A)
log2(max(lambda))
entropy([eps, 1-eps])

N = 1300; 
phi = ( log2 ([1,1] * (A^N) * [1,1]') )/N
phiphi  = log2(max(lambda))

Ido_Ent = entropy([eps, 1-eps]') - log2(max(lambda))
lower_bound = H(p+eps-2*p*eps)
upper_bound = H(p+ 2*eps*(1-eps) -1*p*2*eps*(1-eps))


% Now compute the greatest eignvalue explicitly
lambda = eig(A)
lambda_1 = 0.5* ( (p^2 + (1-p)^2) + sqrt ( (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps*(1-eps)*(eps ^ 2 + (1-eps)^2)) )
lambda_2 = 0.5* ( (p^2 + (1-p)^2) - sqrt ( (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps*(1-eps)*(eps ^ 2 + (1-eps)^2)) )


lamalm = 0.5 * ( A(1,1) + A(2,2) + sqrt( (A(1,1)-A(2,2))^2 + 4*A(1,2)*A(2,1) ) )
A(1,1) + A(2,2) - (p^2 + (1-p)^2)

log2(lambda_1)
