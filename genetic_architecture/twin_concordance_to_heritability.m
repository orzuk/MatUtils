% Compute heritability from twin concordance 
% 
% Input: 
% lambda_mz - mz twins relative risk 
% lambda_dz - dz twins relative risk
% mu - prevalence 
% family_model - how to compute heritability (defaule: ACE)
%
% Output: 
% h_liability - heritabaility on liability scale 
% r_MZ - mz twin concordance on liability scale
% r_DZ - rz twin concordance on liability scale 
% 
function [h_liability r_MZ r_DZ] = ...
    twin_concordance_to_heritability(lambda_mz, lambda_dz, mu, family_model)

if(~exist('family_model', 'var') || isempty(family_model))
    family_model = 'ACE';
end

r_MZ = familial_risk_to_heritability(lambda_mz, 'liability', mu, 1);
%r_DZ = 0.5*familial_risk_to_heritability(lambda_dz, 'liability', mu, 0.5); % Works only for r_DZ <= 0.5
r_DZ = familial_risk_to_heritability(lambda_dz, 'liability', mu, 1);  % just convert to liability scale (no 0.5 here!)

switch family_model % here esitmate heritability from a family by taking DZ and MZ risk
    case 'ACE'
        h_liability = 2*(r_MZ - r_DZ); % standard ACE formula
    case 'ADE'
        h_liability = (4*r_DZ - r_MZ); % /3; % alternative ADE formula
    case 'MZ'
        h_liability = r_MZ;
    case 'DZ'
        h_liability = 2*r_DZ;
    case 'PO'
        h_liability = 2*r_DZ; % assume no dominance 
end


