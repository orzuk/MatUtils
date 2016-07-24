% Compute integral appearing in dalal procedure 
% 
% Input: 
% h - offset 
% k - number of varaibles (look at their maximum) 
% nu - deg. of freedom of t-distrbution 
% 
function ret = two_stage_integral_dalal(h, k, nu) 
    ret = quadgk(@(t) tcdf(t+h, nu).^k .* tpdf(t, nu), -inf, inf); 
    

