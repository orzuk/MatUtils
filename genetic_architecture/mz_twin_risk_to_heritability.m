% Transfer monozigotic twin risk to broad-sense heritability 
% 
% Input: 
% lambda_mz - relative risk for monozigotic twin of affected 
% mu - prevalence of disease
% 
% Output: 
% H - broad sense heritability (on binary scale) 
% 
function H = mz_twin_risk_to_heritability(lambda_mz, mu)

H = mu .* (lambda_mz - 1) ./ (1-mu); 
