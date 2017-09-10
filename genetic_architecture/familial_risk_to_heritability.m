% Convert One familial relative risk to heritability.
% Currently working under the liability-threshold model
%
% Input:
% lambda_R - relative risk for relative of degree R
% h_scale - scale of heritability
% mu - disease prevalence 
% k_R - how close are relatives
%
% Output:
% h - heritability (typically on liability scale) 
% 
function h = familial_risk_to_heritability(lambda_R, h_scale, mu, k_R)

%TOL = 0.00000000000001;
switch h_scale
    case 'binary'
        h = -1; %?? not known yet ..
    case 'liability'
        options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance 
        h = fminbnd(@(x) abs(heritability_to_familial_risk(x, 'liability', mu, k_R)-lambda_R), ...
            0, 1, options); % faster than binary search        
end

