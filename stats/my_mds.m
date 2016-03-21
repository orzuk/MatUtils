% Toy implementation of multidimensional-scaling using eigendecomposition 
% (not efficient and unstable) 
% 
% Input: 
% D - a distance matrix 
% p - dimension 
% 
% Output: 
% V - a set of vectors in R^p with distances approximately D
%
function [V, Q, Lambda] = my_mds( D, p, S)

n = length(D); 
if(~exist('S', 'var') || isempty(S))
    S = repmat(1/n, n, 1); % set default S
    
else
    S = vec2column(S) ./ sum(S); % normalize S 
end
Delta = -0.5 * D.*D; 
C = (eye(n) - repmat(S', n, 1)) * Delta * (eye(n) - repmat(S, 1, n)); 
isposdef_is = isposdef(C)
if(~isposdef(C))
    sprintf('Error !! C Matrix not-positive-definite !!')
end
[Q, Lambda] = eig(C); 
V = Q*sqrt(Lambda); 
V=V(:,1:p); % take top k columns 


