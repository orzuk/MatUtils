% Compute (penalized) log-likelihood for Kronecker product
%
% Input:
% Z - data matrix
% P - first part of covariance matrix
% Q - first part of covariance matrix
% lambda - regularization parameter
%
% Output:
% grad_P - gradient with respect to components of P
%
function [grad_P, grad_Q] = KroneckerLogLikeGradient(Z, P, Q, lambda)

if(~exist('lambda', 'var') || isempty(lambda))
    lambda = 0;
end
p = length(P); q = length(Q); n = size(Z, 1);


grad_P = q*P - Z' * inv(Q) * Z - 4*lambda * inv(P); % New! COmpute Gradient !!
grad_Q = p*Q - Z * inv(P) * Z' - 4*lambda * inv(Q);
