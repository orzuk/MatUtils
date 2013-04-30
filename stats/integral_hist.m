% Perform integral for a histogram (we use simple rectangle integration) 
% 
% Input: 
% x - values
% p - their density
% cumulative_flag - (NEW!) if this is 'on', we return the cumulative function 
% 
% Output: 
% s - the integral  \sum_i p(i) * \delta_x
% p_cum - the cumulative distribution function 
% 
function [s p_cum] = integral_hist(x, p, cumulative_flag)

s = 0.5 .* sum(diff(x) .* (p(1:end-1) +  p(2:end))); 

if(exist('cumulative_flag', 'var') && (~isempty(cumulative_flag)))
    if(cumulative_flag)
        p_cum = [0 0.5 .* (cumsum(diff(x) .* (p(1:end-1) +  p(2:end))))]; 
    end
end