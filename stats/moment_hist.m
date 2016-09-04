% Compute central/non-central moment of any order for a histogram
%
% Input: 
% x - values
% p - their probabilities
% k - the order of moment
% central_flag - compute central moments (default is on) 
% 
% Output: 
% m - moment value: \int_x (x-mu)^k * p(x) dx where mu=0 for central
%     moments and mean for non-central moments 
% 
function m = moment_hist(x, p, k, central_flag, normalize_flag)

if(~exist('central_flag', 'var') || isempty(central_flag))
    central_flag = 1; 
end
if(~exist('normalize_flag', 'var') || isempty(normalize_flag))
    normalize_flag = 1; 
end
if(central_flag) 
    mu = mean_hist(x, p); 
else % compute non-central moment
    mu = 0; 
end
m = sum(p.* (x-mu).^k) / (sum(p)^normalize_flag);
