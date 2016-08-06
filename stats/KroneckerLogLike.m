% Compute (penalized) log-likelihood for Kronecker product
function LL = KroneckerLogLike(Z, P, Q, lambda)

if(~exist('lambda', 'var') || isempty(lambda))
    lambda = 0; 
end
p = length(P); q = length(Q); n = size(Z, 2); 

Sigma = kron(P, Q); 
LL =  (n/2) * ( logdet(P) + logdet(Q) - trace(S/Sigma) ) + ... % Add penalized part
    lambda * ( norm(inv(P), 'fro') + norm(inv(P), 'fro') ); 


           