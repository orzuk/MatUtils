% Given p-values and R (num rejections) compute q
% used to choose the cutoff such that exactly R will be rejected  
%
% Input: 
%  P      matrix of SORTED(!) p-values (each row is a vector)
%  R      vector containing indices of hypotheses that passsed selected procedure
%  str    name of procedure (case insensetive)
%         function also outputs V - the number of false positives. Note: We
%         assume that thet first m0 pvals are the null
%  is_sorted  flag saying that the p-values are already sorted, to save  sorting time (default is 'false')
%  TOL    at what tolerance do we want to find q   
%  Output:
%  q      scalar controlling FDR output
%
function [q q_s] = q_from_R_FDR_mat_main(P, R, str, is_sorted, TOL, varargin)
if(~exist('TOL', 'var')) % set default tolerance
    TOL = 0.0000001;
end

if(~exist('is_sorted', 'var'))
   is_sorted = 0; 
end
if(~is_sorted) % sort for first time
    P = sort(P,2); 
end

% Naive way: Perform bisection method (binary search) and get the maximal q which rejects R
q_min = 0; q_max = 1; 
while((q_max - q_min) > TOL)
    q_mid = (q_min+q_max)/2;
    R_mid = FDR_mat_main(P, q_mid, str, [], [], 1); % assume sorted
    if(R_mid < R)
        q_min = q_mid;
    else
        q_max = q_mid;
    end
end
q = q_mid; % return the mid value 

[n,m] = size(P); 
q_s = P .* m ./ repmat([1:m], n, 1);
for i=n-1:-1:1
    q_s(i,:) = min(q(i,:), q_s(i+1,:));
end
q_s = q_s(:,max(R,1)); q_s(R == 0) = 0; % return the R's q-value


