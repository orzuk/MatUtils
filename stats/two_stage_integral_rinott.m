% Compute integral appearing in dalal procedure 
% 
% Input: 
% h - offset 
% k - number of varaibles (look at their maximum) 
% nu - deg. of freedom of t-distrbution 
% 
% Output: 
% ret - [int_{t=0}^{infinity} G_{nu}(t+h) * g_{nu}(t) dt]^k 
% where G_{nu}, g_{nu} are t-student c.d.f. and p.d.f.
% 
function ret = two_stage_integral_rinott(h, k, nu) 
    ret = quadgk(@(t) tcdf(t+h, nu) .* tpdf(t, nu), -inf, inf).^k; 
    

