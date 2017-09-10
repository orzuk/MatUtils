% Compute a 'power-index' saying how powerful a given sample is. 
% This can be used to compare case-control studies with different
% proportions, as well as QTL's. 
% This number is NOT power, but non-centrality parameter indicative of power
% 
% Input: 
% num_cases - # of cases (binary trait) or samples (QTL)
% num_controls - # of controls (binary trait)
% trait_type - binary or QTL
% mu - disease prevalence (for binary traits) 
% 
% Output: 
% power_effect_index - index represening non-centrality power parameter
% 
function power_effect_index = sample_size_to_power_effect_index(num_cases, num_controls, ...
    trait_type, mu)

f_vec = 0.5; h_liab = 0.0000001;  % set arbitrary effect size and RAF (explain 1% of variance)
if(~exist('mu', 'var') || isempty(mu))
    mu = 0.3; % prevalence in case-control study 
end
switch trait_type
    case {'binary', 'Binary'}
        grr_vec = heritability_to_genetic_relative_risk(h_liab, 'liability', f_vec, mu);
        power_effect_index = genetic_relative_risk_to_non_centrality_parameter( ... 
            f_vec, grr_vec, num_cases, num_controls, mu, 'ott'); 
    case {'QTL', 'Quantitative'}
        power_effect_index = variance_explained_to_non_centrality_parameter(h_liab, num_cases); 
end
power_effect_index = power_effect_index * (1-h_liab)/h_liab;
        



