% Perform cumulative summation for a histogram 
% (we use simple rectangle integration. Last value is integral_sum) 
% 
% Input: 
% x - values
% p - their density
% 
% Output: 
% c - the cumsumintegral  \sum_i p(i) * \delta_x
%
function c = cumsum_hist(x, p)

%c = 0.5 .* cumsum(diff(x) .* (p(1:end-1) +  p(2:end))); 

d = diff(x); 

c = cumsum(p .* 0.5 .* ([d 0] + [0 d])); % give edges smaller weight (half) 


