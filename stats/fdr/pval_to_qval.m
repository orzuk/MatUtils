% Computes for a given p-vals vector and a given R a bound q on the FDR.
% The method is by computing
% first the standard q-value (from the standard BH procedure), and then
% correct it using the procedure-specific m0 estimator
%
% Input: 
%   P - Vector of matrix of SORTED(!) p-values (each row is a vector)
%   (optional) is_sorted - flag saying if the p-values are already sorted, to save sorting time - Default is 'false'
%
% Output:
%   q - the q-values
%
function q = pval_to_qval(P, is_sorted, varargin)

if(~exist('is_sorted', 'var'))
   is_sorted = 0; 
end
if(~is_sorted) % sort for first time
    [P, sort_perm] = sort(P,2); 
end

[n,m] = size(P); 
q = min(P .* m ./ repmat(1:m, n, 1), 1);
for i=m-1:-1:1
    q(:,i) = min(q(:,i), q(:,i+1));
end
for j=1:n
    q(j,:) = q(j,inv_perm(sort_perm(j,:)));
end



