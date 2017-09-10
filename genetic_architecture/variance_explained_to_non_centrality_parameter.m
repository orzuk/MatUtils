% Compute chi-square non-centrality parameter for a QTL locus (indicative of power)
% Problem: this depends on the effect size (fraction of variance explained). There's no 'universal' power measure
% 
% Input: 
% V - effect size (fraction of phenotypic variance explained). This indirectly depends on the risk-allele-frequency
% num_samples - # of cases in study
% 
% Output: 
% non_centrality_parameter - chi-square non-centrality. Marks deviation from null-hypothesis and indicates power
%
function non_centrality_parameter = variance_explained_to_non_centrality_parameter( ...
    V, num_samples)

non_centrality_parameter = num_samples .* V ./ (1-V); %     NCP <- N*q.sq/(1-q.sq);
