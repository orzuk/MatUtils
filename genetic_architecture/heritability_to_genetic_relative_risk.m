% Compute relative risk from heritability (either on binary or liability scale)
%
% Input:
% h - heritability value
% h_scale - what scale was heritability measured on (e.g. 'binary' vs 'liability')
% f_vec - minor allele frequencies
% mu - disease prevalence
%
% Output:
% GRR - genetic relative risk for each locus (computed using an approximation)
% GRR_exact - same GRR, computed exactly using Gaussian integrals
% 
function [GRR GRR_exact] = heritability_to_genetic_relative_risk(h, h_scale, f_vec, mu)

switch h_scale
    case 'liability' % Perform an additional transformation
        h_liab=h; % save original heritability 
        h = heritability_scale_change(h, 'binary', mu); % [h, ~, h_exact]  = .. convert to binary. This uses Falconer's approximation    
    case 'binary' % here do nothing
        
end
%h, h_exact, h-h_exact

tmp_const = f_vec.^2 - mu.*f_vec.*(1-f_vec) ./ (h.*(1-mu)); % an expression appearing twice in the GRR formula 
GRR = (-f_vec - sqrt(f_vec.^2 - tmp_const)) ./ tmp_const + 1; % Copmute approximation. We took negative values so subtract sqrt 

if(nargout>1) % New! compute exact GRR using Gaussian integrals
    h = h_liab;
    x_mu = norminv(1-mu);    
    beta = sqrt(h ./ (f_vec.*(1-f_vec)));
    GRR_exact = ( 1-normcdf( (x_mu-beta.*(1-f_vec)) ./ sqrt(1-h) ) ) ./ ...
        ( 1-normcdf( (x_mu+beta.*f_vec) ./ sqrt(1-h) ) );    
end

