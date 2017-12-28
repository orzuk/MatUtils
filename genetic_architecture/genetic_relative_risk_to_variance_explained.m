% Compute heritability on liability scale assuming an additive/multiplicative (?) model
%
% Input:
% f_vec - risk allele frequencies
% grr_vec - genetic relative risk for each SNP
% mu - disease prevalence
% SNP_TYPE - diploid (default) or binary
% model - (optional) multiplicative or additive (or logistic?)
% scale (new) report variance on what scale 
% 
% Output:
% h_liab - heritability on liability scale
% h_liab_vec - heritability on liability scale of each locus 
% 
function [h_liab, h_liab_vec] = ...
    genetic_relative_risk_to_variance_explained(f_vec, grr_vec, mu, SNP_TYPE, scale)

if(~exist('SNP_TYPE', 'var') || isempty(SNP_TYPE))
    SNP_TYPE = 'diploid';
end
if(~exist('scale', 'var') || isempty(scale))
    scale = 'liability';
end
AssignGeneralConstants;
f_vec = vec2column(f_vec);

if(length(grr_vec) > length(f_vec))
    f_vec = repmat(f_vec, 1, length(grr_vec));
end

switch SNP_TYPE % Find impossible inds
    case {'binary', 'haploid'}
        bad_inds = find( (grr_vec > (1-f_vec) ./ (mu-f_vec)) & (f_vec < mu) );
    case 'diploid' % more complicated equation: (f^2-mu)*GRR^2 + 2*f*(1-f)*GRR+(1-f)^2 >= 0
        bad_inds = find( (f_vec.^2 - mu).*grr_vec.^2 + 2.*f_vec.*(1-f_vec).*grr_vec + (1-f_vec).^2 < 0 ); 
%         bad_inds = find( (grr_vec < (sqrt(mu)-f).*(1-f_vec) ./ (f_vec.^2-mu.^2)) & (f_vec.^2 < mu) );    % 
%         bad_inds = union(bad_inds, ...
%             find( (grr_vec > (sqrt(mu)+f).*(1-f_vec) ./ (f_vec.^2-mu.^2)) & (f_vec.^2 < mu) ) ); 
%         bad_inds = union(bad_inds, ...
%             find( 9999*(grr_vec > (sqrt(mu)+f).*(1-f_vec) ./ (f_vec.^2-mu.^2)) & (f_vec.^2 > mu) ) ); % here a < 0 
        
        
        
        % find( (grr_vec > (1-f_vec) .* ( -f_vec + sqrt(f_vec.^2 - (4-mu)^2) ) ./ (1-mu)) & (f_vec.^2 < mu) );         
end

V_mult = mu.*(1-mu); % for additive heritability we need also prevalence
switch SNP_TYPE
    case {'binary', 'haploid'}
        beta_vec = mu .* (grr_vec - 1) ./ (1 - f_vec + f_vec .* grr_vec); % beta_vec = grr_vec;
        V_add = f_vec .* (1-f_vec) .* beta_vec.^2; % new: marginal additive variance in additive model
    case 'diploid'
        beta_vec = mu ./ ...
            ( (1-f_vec).^2 + 2.*f_vec.*(1-f_vec) .* grr_vec + f_vec.^2 .* grr_vec.^2 );
        V_add = f_vec .* (1-f_vec) .* beta_vec.^2 .* ...
            ( (1-f_vec).*(2-f_vec) - 4.*(1-f_vec).^2.*grr_vec + 2.*(1-3.*f_vec.*(1-f_vec)).*grr_vec.^2 - ...
            4.*f_vec.^2.*grr_vec.^3 + f_vec.*(1+f_vec).*grr_vec.^4 ); % new: marginal additive variance in additive model with two alleles 

%         % Alternative: compute the actual probability distribution
%         W = (1-f_vec).^2 + 2.*f_vec.*(1-f_vec) .* grr_vec + f_vec.^2 .* grr_vec.^2;
%         p_x_z = zeros(2,3); 
%         p_x_z = repmat([ (1-f_vec).^2 2.*f_vec.*(1-f_vec) f_vec.^2 ], 2, 1) .* ...
%             [(mu .* grr_vec .^ [0 1 2] ./ W)' (1 - mu .* grr_vec .^ [0 1 2] ./ W)']'
%         V_add3 = sum( p_x_z(1,:) .* p_x_z(2,:)  ./ sum(p_x_z) );
%         V_add3 = V_mult - V_add3
%         
%         V_add4 = mu^2 .* [ (1-f_vec).^2 + 2.*f_vec.*(1-f_vec) .* grr_vec.^2 + f_vec.^2 .* grr_vec.^4  - W^2] ./ W.^2

end
h_add = sum(V_add) ./ V_mult; % if we're given multiple mu's this will give a vector !!!



switch scale
    case 'liability' % default
        h_liab = heritability_scale_change(h_add, 'liability', mu); %h_add * V_mult / z_score^2; % Compute also liability threshold model heritability
        h_liab_vec = heritability_scale_change(V_add ./ V_mult, 'liability', mu); % liability-scale heritability of each variant
    case 'binary'
        h_liab = h_add;
        h_liab_vec = V_add ./ V_mult;
end
% switch SNP_TYPE % just multiply by two variance explained 
%     case 'diploid'
%         h_liab = h_liab .* 2; h_liab_vec = h_liab_vec .* 2;
% end

if(~isempty(bad_inds))
    sprintf('Error! GRR cannot be that large!!!')
    h_liab_vec(bad_inds) = NaN; h_liab = NaN;
end
