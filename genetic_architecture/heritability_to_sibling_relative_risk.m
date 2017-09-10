% Compute sibling and mz twins relative risks from heritability
% 
% Input: 
% h - heritability
% h_scale - scale of heritability 
% mu - trait mean (prevalence) 
%
% Output: 
% lambda_s - sibling relative risk
% lambda_mz - monozigotic twin relative risk
% 
function [lambda_s lambda_mz] = ... %  lambda_s_mult lambda_s_general lambda_mz_general] = ...
    heritability_to_sibling_relative_risk(h, h_scale, mu)
TOL = 0.00000000000001;
threshold = norminv(1-mu); z_score = normpdf(threshold); % Get Gaussian height at the incidence threshold
num_inputs = max(length(threshold), length(h)); 
if(length(threshold) == 1)
    threshold = repmat(threshold, 1, num_inputs); 
end
if(length(mu) == 1)
    mu = repmat(mu, 1, num_inputs); 
end
if(length(h) == 1)
    h = repmat(h, 1, num_inputs); 
end

switch h_scale
    case 'binary'
        H = h; % don't change heritability
    case 'liability'
        [~, H] = heritability_scale_change(h, 'binary', mu);
end
lambda_mz = heritability_to_mz_twin_risk(H, mu);
lambda_s = zeros(1,length(h));
switch h_scale
    case 'binary'
        lambda_s = sqrt(lambda_mz); % a (possibly bad approximation)
    case 'liability'
        for i=1:length(h)
            lambda_s(i) = quadl(@(x)heritability_to_familial_risk_integrand(x, threshold(i), sqrt(h(i)), 0.5), ...
                threshold(i), threshold(i)+10*max(h(i),1-h(i)), TOL) ./ mu(i)^2; % compute EXACT sibling relative risk under the liability threshold model
        end
end

% lambda_s = lambda_s_mult; % quadl(@(x)heritability_to_sib_fun(x, threshold, h), -10*max(h,1-h), 10*max(h,1-h), TOL) ./ mu^2; % compute EXACT sibling relative risk under the liability threshold model
% lambda_mz_general = quadl(@(x)heritability_to_familial_risk(x, threshold, sqrt(h), 1), ...
%     threshold, threshold+10*max(h,1-h), TOL) ./ mu^2; % compute EXACT mz twin relative risk under the liability threshold model
% lambda_mz3 = quadl(@(x)heritability_to_lambda_mz_local(x, threshold, sqrt(h)), ...
%     -threshold-10*max(h,1-h), threshold+10*max(h,1-h), TOL) ./ mu^2;


