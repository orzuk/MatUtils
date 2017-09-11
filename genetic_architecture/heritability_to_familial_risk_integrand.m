% Internal function computing integrand for converting heritability to familial risk based on joint gaussian distribution
% 
% Input: 
% x - value on liability scale 
% T - liability threshold 
% h - heritability on liability scale 
% k_R - relative IBD sharing of relatives (this is twice the kinship coefficient)
% 
% Output: 
% y - value used for the familial risk integral formula 
% 
function y = heritability_to_familial_risk_integrand(x, T, h, k_R)

tmp_sigma = sqrt(1 - h.^4.* k_R.^2);
y = normpdf(x) .* (1 - normcdf( ((T - h.^2.*k_R.*x) ./ tmp_sigma) ));

