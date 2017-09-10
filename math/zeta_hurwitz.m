% Compute the Hurwitz Zeta function
% This is a quick-dirty way to compute the function by 
% trancating a power-series. It assumes that the series converges.
% Written since Matlab currently doesn't have an implementation (only for
% the standard Rieman Zeta function)
%
% Input: 
% s - first argument (must have Re(s) > 1)
% q - parameter (must have Re(q) > 0) 
% max_n - (optional) where to truncate the power series
% 
% Output: 
% z - value of the Hurwitz Zeta function
% 
function z = zeta_hurwitz(s, q, max_n)

num_vals = max(length(s), length(q)); z = zeros(num_vals,1);
if(length(s) == 1)
    s = repmat(s,num_vals,1);
end
if(length(q) == 1)
    q = repmat(q,num_vals,1);
end

% % % % if(~exist('max_n', 'var') || isempty(max_n)) % set maximum n to sum (determines accuracy) 
% % % %     max_n = 10000; 
% % % % end
% % % % 
% % % % for i=1:num_vals
% % % %     z(i) = sum(1 ./ (q(i) + 0:max_n).^s(i)); % This works when Re(s) > 1. What happens when it's smaller? 
% % % % end

% Alternative: integral representation
for i=1:num_vals
    z(i) = 1 ./ (2.*q(i).^s(i)) + q(i).^(1-s(i)) ./ (s(i)-1) + 2 .* quadl(@(x)internal_zeta(x, s(i), q(i)), 0,  9999);
end

function ret = internal_zeta(x, s, q) 

ret = sin(s .* atan(x ./ q)) ./ ( (exp(2*pi.*x)-1) .* (q.^2+x.^2).^(s./2) );



