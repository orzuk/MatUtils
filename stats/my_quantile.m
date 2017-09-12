% Like Matlab's quantile but enables to get quantiles outside [0,1].
% So one can for example ask for the -5% quantile or 150% quantile
% 
% Input: 
% V - vector of values
% alpha - the fraction we want
% 
% Output: 
% q - the quantile at fraction alpha
function q = my_quantile(V, alpha)
if( (alpha >= 0) && (alpha <= 1) )
    q = quantile(V, alpha);
else
    alpha_frac = alpha - floor(alpha);
    alpha_div = floor(alpha);
    V_range = max(V) - min(V);
    q = quantile(V, alpha_frac) + V_range *  alpha_div;
end

