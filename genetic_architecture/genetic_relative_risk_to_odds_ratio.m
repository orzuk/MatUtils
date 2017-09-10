% Convert genetic relative risk to odds-ratio
% 
% Input: 
% GRR - 
% f_vec - allele frequency
% mu - prevalence 
% 
function OR = genetic_relative_risk_to_odds_ratio(GRR, f_vec, mu)

% OR = (1 - 2.*f_vec - GRR .* f_vec + f_vec.^2 + GRR.*f_vec.^2 - mu + mu.*f_vec) ./ ...
OR = GRR .* (1 - f_vec + GRR .* f_vec - mu) ./ ...
(1 - f_vec + GRR .* f_vec - GRR .* mu);


