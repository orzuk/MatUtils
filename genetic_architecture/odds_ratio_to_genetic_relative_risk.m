% Convert odds-ratio to genetic relative risk
% 
% Input: 
% OR - odds ratio for each locus
% CAF - Control Allele Frequency for each locus
% mu - disease prevalence
% 
% Output: 
% GRR - disease relative risk for each locus
% RAF - risk allele frequency for each locus 
% 
function [GRR, RAF] = odds_ratio_to_genetic_relative_risk(OR, CAF, mu)

%GRR = (1 - 2.*f_vec + f_vec.^2 - mu + mu.*f_vec + OR .*(f_vec-1)) ./ ...
%    (f_vec - f_vec.^2 - OR .* (f_vec + mu));


RAF = CAF .* (1 - mu + mu .* CAF ./ (1 - CAF + OR .* CAF)); % set risk-allele frequency in the population! 

GRR = ( OR .* RAF - OR .* mu + mu + RAF - 1 + ...
    sqrt( (1-RAF - mu + OR .* mu - OR.*RAF).^2 - 4.*OR.*RAF.*(RAF-1) ) ) ./ (2 .* RAF);



