% Computes for a given p-vals vector and a given R a bound q on the FDR.
% The method is by computing
% first the standard q-value (from the standard BH procedure), and then
% correct it using the procedure-specific m0 estimator
%
% Input: 
%   P - Vector of matrix of SORTED(!) p-values (each row is a vector)
%   R - Number of hypotheses rejected
%   
%   (optional) proc_str -  name of procedure (case insensetive).
%       Default procedure: 
%           'ibh_up' (Improved Benjamini Hochberg)
%       Other Supported procedures: 
%           'bh95' (Standard Benhamini Hochberg)
%           'sth' (Storey's procedure with lambda=0.5)
%   (optional) is_sorted - flag saying if the p-values are already sorted, to save sorting time - Default is 'false'
%
% Output:
%   q - The estimated FDR of the R lowest p-values
%
%   Example: 
%   P = rand(1,1000); % generate a vector of p-values, of which a 100 are taken from the alternative hypothesis
%   P(1:100) = P(1:100) ./ 200; 
%   q = FDR_R_to_q(P, 100)  % Compute the FDR of the lowest 100 p-values
%
function q = FDR_R_to_q(P, R, proc_str, is_sorted, varargin)

if(~exist('proc_str', 'var'))
    proc_str = 'ibh_up';
end
if(isempty(proc_str))
    proc_str = 'ibh_up';
end
if(~exist('is_sorted', 'var'))
   is_sorted = 0; 
end
if(~is_sorted) % sort for first time
    P = sort(P,2); 
end

[n,m] = size(P); 
q = min(P .* m ./ repmat(1:m, n, 1), 1);
for i=n-1:-1:1
    q(i,:) = min(q(i,:), q(i+1,:));
end
q = q(:,max(R,1)); q(R == 0) = 0; % return the R's q-value (BH95)

q = q .* m0_estimator(P, proc_str) ./ m; % Adjust the q to the desired procedure


