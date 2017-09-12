% Compute integral appearing in dalal procedure
%
% Input:
% h - offset
% k - number of varaibles (look at their maximum)
% nu - deg. of freedom of t-distrbution
% limit - where to take integration (default from -infinity to infinity)
%
% Output:
% ret - [int_{t=0}^{infinity} G_{nu}(t+h) * g_{nu}(t) dt]^k
% where G_{nu}, g_{nu} are t-student c.d.f. and p.d.f.
%
function ret = two_stage_integral_rinott(h, nu, limit)

if(~exist('limit', 'var') || isempty(limit))
    limit = [-inf inf];
end
if(isscalar(limit))
    limit = [-limit limit];
end

if(nu == Inf) % Gaussian
    %%    ret = quadgk(@(t) normcdf(t+h) .* normpdf(t), limit(1), limit(2)) .^k; % change boundaries from -infinity to infinity
    ret = quadgk(@(t) (normcdf(-t-h)) .* normpdf(t), limit(1), limit(2)); % change boundaries from -infinity to infinity
else
    %%    ret = quadgk(@(t) (tcdf(-t-h, nu)) .* tpdf(t, nu), limit(1), limit(2)); %  .^k; % change boundaries from -infinity to infinity
    ret = quadgk(@(t) (tcdf(-t-h, nu)) .* tpdf(t, nu), limit(1), 0, 'AbsTol', 10^(-20), 'RelTol', 10^(-13)) + ... %  .^k; % change boundaries: assume first is negative and second positive !!
        quadgk(@(t) (tcdf(-t-h, nu)) .* tpdf(t, nu), 0, limit(2), 'AbsTol', 10^(-20), 'RelTol', 10^(-13)); %  .^k; %
end

