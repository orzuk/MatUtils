% Transfer broad-sense heritability to monozigotic twin risk
% 
% Input: 
% H - broad sense heritability
% mu - prevalence
% 
% Output: 
% lambda_mz - relative risk for twin
% 
function lambda_mz = heritability_to_mz_twin_risk(H, mu)

lambda_mz = 1 + H .* (1-mu) ./ mu;
