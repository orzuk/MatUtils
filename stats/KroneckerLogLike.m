% Compute (penalized) log-likelihood for Kronecker product
%
% Input:
% Z - data matrix
% P - first part of covariance matrix
% Q - first part of covariance matrix
% lambda - regularization parameter
%
function LL = KroneckerLogLike(Z, P, Q, lambda)

if(~exist('lambda', 'var') || isempty(lambda))
    lambda = 0;
end
p = length(P); q = length(Q); n = size(Z, 1);
S = (Z'*Z)./ n; % Get sample covariance matrix

Sigma = kron(P, Q);
if(isscalar(lambda))
    lambda = [lambda lambda];
end
LL =  -(n/2) * ( q*logdet(P) + p*logdet(Q) + trace(S/Sigma) ) -  ... % Add penalized part
    lambda(1) * norm(inv(P), 'fro')^2 + lambda(2) * norm(inv(Q), 'fro')^2;


% LL = -lambda * ( norm(inv(P), 'fro')^2 ); % -(n/2) * ( q*logdet(P) ) %  -(n/2) * (  trace(S/Sigma) )

