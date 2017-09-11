% Convert sibling-relative risk to heritability 
% under the assumption of a multiplicative model (is this valid?).
function h = sibling_relative_risk_to_heritability(lambda_s, h_scale, mu)

lambda_mz = lambda_s.^2; % take monozigotic-twin risk 
H = mz_twin_risk_to_heritability(lambda_mz, mu); % broad sense heritability

switch h_scale
    case 'binary' % don't change heritability 
        h = H; 
    case 'liability' % move to liability scale
        [~, h] = heritability_scale_change(H, 'liability', mu);
end




