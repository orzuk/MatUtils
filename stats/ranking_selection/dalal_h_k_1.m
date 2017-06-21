% Compute asymptotic of hk1
%
% Input:
% nu - degrees of freedom (N0-1 where N0 is initial sample size)
% k - # of populations
% p - desired pcs
% derivative: 0 - compute h, 1 - compute first derivative of log, 2 -
% compute 2nd derivative without (1/2*nu^2) part
%
% Output:
% h - result
%
function h = dalal_h_k_1(nu, k, p, der_m)

q_p = (-1./log(p)).^(1./nu);

if(~exist('der_m', 'var') || isempty(der_m))
    der_m = 0;
end
switch der_m
    case 0
%        h = (gamma((nu+1)/2) / (nu^(1-nu/2)*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p; % take asymptotic solution
        h = exp( (1./nu) .* (  log( -k ./ (sqrt(pi) .* log(p)) ) + gammaln((nu+1)./2) - gammaln(nu./2) + (nu./2-1).*log(nu) ) ); 
    case 1 % take 1st derivatove of log(h)
        h = (1./(2.*nu.^2)) .* ( -2 - 2*log( -k ./ (sqrt(pi)*log(p)) ) + nu + 2*log(nu) + 2 * (log (gamma(nu./2)) - log (gamma((nu+1)./2))) ...
            + (nu./2) .* (psi(1,(nu+1)./2)-psi(1,nu./2)) );
    case 2 % take derivative of H
        h = 1 + 2./nu - 0.5.*(psi(1,nu./2)-psi(1,(nu+1)./2));
end

syms x; 

