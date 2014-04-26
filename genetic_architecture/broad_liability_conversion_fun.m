
% Function is: (phi(x/h)/h) * (1 - Phi((T-x)/sqrt(1-h^2)^2)
function y = broad_liability_conversion_fun(x, T, h) % compute conditional variance of Z on X for broad sense heritability

% h = sqrt(h); % take h and not heritability (already taken care of in input!)
y = ( exp(-x.^2/(2.*h.^2)) ./ (h.*sqrt(2.*pi)) ) .* (1-normcdf((T-x) ./ sqrt(1-h.^2))).^2; % 1-h.^2

