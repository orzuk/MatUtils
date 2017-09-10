% Convert One familial relative risk to heritability.
% Currently working under the liability-threshold model
%
% Input:
% lambda_R - relative risk for relative of degree R
% h_scale - scale of heritability
% h - heritability (typically on liability scale) 
% k_R - how close are relatives
%
% Output:
% mu - disease prevalence 
% 
function mu = familial_risk_and_heritability_to_prevalence(lambda_R, h_scale, h, k_R)

switch h_scale
    case 'binary'
        mu = -1; %?? not known yet ..
    case 'liability'
        
        
        mu = fminbnd(@(x) abs(heritability_to_familial_risk(h, 'liability', x, k_R)-lambda_R), 0, 1); % faster than binary search        
%         h_min = 0; h_max = 1; h_mid = 0.5;
%         while (h_max - h_min > TOL)
%             lambda_R_output = heritability_to_familial_risk(h_mid, 'liability', mu, k_R); 
%             if(lambda_R_output > lambda_R)
%                 h_max = h_mid;
%             else
%                 h_min = h_mid;
%             end
%             h_mid = (h_max + h_min)/2; %  h_binary_broad_output
%         end
%         h = h_mid;        
end

