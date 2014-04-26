% Compute relative risks for a family member from heritability 
% under LT model
%
% Input:
% h - heritability
% h_scale - scale of (input) heritability
% mu - trait mean (prevalence)
% k_R - amound of IBD shared (this is twice the kinship coefficient: 1 for MZ, 0.5 for DZ or sibs and so on)
%
% Output:
% lambda_R - relative risk for a family member of degree R
%
function lambda_R  = heritability_to_familial_risk(h, h_scale, mu, k_R)

TOL = 0.00000000000001;
threshold = norminv(1-mu); % z_score = normpdf(threshold); % Get Gaussian height at the incidence threshold
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

% We don't really need big H !!!
% switch h_scale
%     case 'binary'
%         H = h; % don't change heritability
%     case 'liability'
%         [~, H] = heritability_scale_change(h, 'binary', mu);
% end

lambda_R = zeros(1,length(h));
switch h_scale
    case 'binary'
        lambda_R = -1; % NOT WORKING YET! sqrt(lambda_mz); % a (possibly bad approximation)
    case 'liability'
        for i=1:length(h)
            if(h(i)==0) % special case (no risk - singular case) 
                lambda_R(i)=1;
            else
                lambda_R(i) = quadl(@(x)heritability_to_familial_risk_integrand(x, threshold(i), sqrt(h(i)), k_R), ...
                    threshold(i), threshold(i)+10*max(h(i),1-h(i)), TOL) ./ mu(i)^2; % compute EXACT familial relative risk under the liability threshold model
            end
        end
end

