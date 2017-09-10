% Compute integral appearing in dalal procedure 
% 
% Input: 
% h - offset 
% k - number of variables (look at their maximum) 
% nu - deg. of freedom of t-distrbution 
% 
% Output: 
% ret - int_{t=0}^{infinity} G_{nu}(t+h)^k * g_{nu}(t) dt 
% where G_{nu}, g_{nu} are t-student c.d.f. and p.d.f.
% 
function ret = two_stage_integral_dalal(h, k, nu)

if(isscalar(nu)) % here we have k1, k2
    nu = [nu nu]; 
end
    
    
if(~isscalar(k)) % here we have k1, k2
    ret = quadgk(@(t) tcdf(h-t, nu(1)).^k(1) .* k(2) .* tcdf(t, nu(2)).^(k(2)-1) .* tpdf(t, nu(2)), -inf, inf, 'AbsTol', 10^(-15), 'RelTol', 10^(-10));
else
    if(nu(1) == Inf) % take Gaussian approximation
        ret = quadgk(@(t) normcdf(t+h).^k .* normpdf(t), -inf, inf, 'AbsTol', 10^(-15), 'RelTol', 10^(-10));
    else
        ret = quadgk(@(t) tcdf(t+h, nu(1)).^k .* tpdf(t, nu(2)), -inf, inf, 'AbsTol', 10^(-15), 'RelTol', 10^(-10));
    end
end

% New: take derivative of g(nu) to find best nu: 
%(((c x^(-1 + x/2) Gamma[(1 + x)/2])/Gamma[x/2])^x^(-1) (-2 + x + x Log[x] - 
%2 Log[(c x^(-1 + x/2) Gamma[(1 + x)/2])/Gamma[x/2]] - x PolyGamma[0, x/2] + x PolyGamma[0, (1 + x)/2]))/(2 x^2)

% take derivative of log(g)