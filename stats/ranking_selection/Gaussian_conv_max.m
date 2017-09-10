% Compute integral appearing in dalal procedure 
% 
% Input: 
% h - offset 
% k - number of variables (look at their maximum) 
% nu - deg. of freedom of t-distrbution 
% 
% Output: 
% ret - int_{t=-infinity}^{infinity} s * Phi(t+h)^(k-s) * Phi(t)^(s-1) * phi(t) dt 
% where Phi, phi are standard Gaussian c.d.f. and p.d.f.
% 
% 
% Compute PCS: Prob( max_{i=1..k-s} X_i <= min (j=1..s) Y_j - Delta ) for
% X_i, Y_j ~ N(0,1) i.i.d. 
function ret = Gaussian_conv_max(Delta, k, s)

% take Gaussian integral 
ret = quadgk(@(t) normcdf(-t+Delta).^(k-s) .* s .* normcdf(t) .^ (s-1) .* normpdf(t), -inf, 0, 'AbsTol', 10^(-17), 'RelTol', 10^(-11)) + ...
    quadgk(@(t) normcdf(-t+Delta).^(k-s) .* s .* normcdf(t) .^ (s-1) .* normpdf(t), 0, inf, 'AbsTol', 10^(-17), 'RelTol', 10^(-11));


