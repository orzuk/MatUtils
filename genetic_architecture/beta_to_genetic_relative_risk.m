% Convert beta on liability scale to grr (on disease scale)
% 
% Input: 
% beta_vec - effect size (regression coefficient)
% f_vec - allele freq. 
% mu - prevalence
%
% Output: 
% grr_vec - genetic relative risk
% grr_vec_binary - genetic relative risk when the beta is associated with a binary random variable (not Gaussian)
%
function [grr_vec grr_vec_binary risk_given_risk_allele risk_given_harmless_allele prevalence_vec] = ...
    beta_to_genetic_relative_risk(beta_vec, f_vec, mu)

num_beta = length(beta_vec);
if(length(f_vec) == 1)
    f_vec = repmat(f_vec, size(beta_vec));
end
neg_inds = find(beta_vec < 0); 
beta_vec(neg_inds) = -beta_vec(neg_inds); 
f_vec(neg_inds) = 1-f_vec(neg_inds);
h_vec = 2 .* f_vec .* (1-f_vec) .* beta_vec.^2;
good_inds = 1:num_beta;
if(max(h_vec) > 1) 
    bad_inds = find(h_vec > 1);
    good_inds = setdiff(1:num_beta, bad_inds); 
%    error('Error! beta and f values lead to heritability > 100% !!!'); 
    sprintf('Error! beta and f values lead to heritability > 100% !!!') 
else
    bad_inds = []; 
end
[grr_vec(good_inds) grr_vec_exact(good_inds)] = heritability_to_genetic_relative_risk(h_vec(good_inds), 'liability', f_vec(good_inds), mu); 
if(max(h_vec) > 1) 
    grr_vec(bad_inds) = NaN; % 99999999999999999;
    grr_vec_exact(bad_inds) = NaN; % 99999999999999999;    
end
grr_vec(neg_inds) = 1 ./ grr_vec(neg_inds); % flip grr for negative beta 



% New: compute directly (faster and should give correct results - good for
% debugging)
x_mu = norminv(1-mu);
%beta_vec = sqrt(2) .* beta_vec; % make everything diploid
denom = sqrt(1 - 2.*beta_vec.^2 .* f_vec .* (1-f_vec)); % smaller variance in the denominator 
risk_given_risk_allele = ( 1 - normcdf(   (x_mu - sqrt(2) .* beta_vec(good_inds) .* (1-f_vec(good_inds))) ./ denom(good_inds) ) );
risk_given_harmless_allele = ( 1 - normcdf(   (x_mu + sqrt(2) .* beta_vec(good_inds) .* f_vec(good_inds)) ./ denom(good_inds) ) );

grr_vec_binary(good_inds) = risk_given_risk_allele ./ risk_given_harmless_allele;
% ( 1 - normcdf(   (x_mu - sqrt(2) .* beta_vec(good_inds) .* (1-f_vec(good_inds))) ./ denom(good_inds) ) ) ./ ...
%     ( 1 - normcdf(   (x_mu + sqrt(2) .* beta_vec(good_inds) .* f_vec(good_inds)) ./ denom(good_inds) ) );
if(max(h_vec) > 1) 
    grr_vec_binary(bad_inds) = NaN; % 99999999999999999;
end


prevalence_vec(good_inds) = risk_given_risk_allele .* f_vec(good_inds) + ...
    risk_given_harmless_allele .* (1-f_vec(good_inds));
prevalence_vec(bad_inds) = NaN;
