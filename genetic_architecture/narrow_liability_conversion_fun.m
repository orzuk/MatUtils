% Compute liability function : x*phi(x) * (1 - Phi((T-h_x*x)/sqrt(1-h^2))
% 
% Input: 
% x - value of liability 
% x_mu - threshold 
% h - how much heritability by liability 
%
% Output: 
% y - density for computing E[xz]
%
function y = narrow_liability_conversion_fun(x, x_mu, h_x) % compute correlation between Z and X for narrow sense heritability

try_func = 1;
switch try_func
    case 0 % OLD function - probably wrong!!!
        y = h_x .* ( x.*exp(-x.^2/(2)) ./ (sqrt(2.*pi)) ) .* (1-normcdf((x_mu-x) ./ sqrt(1-h.^2))); % 1-h.^2
    case 1 % NEW function - Working !!!! 
        y = ( x.*exp(-x.^2/2) ./ (sqrt(2.*pi)) ) .* ...
            normcdf(-(x_mu-h_x.*x) ./ sqrt(1-h_x.^2)); % 1-h.^2
end
