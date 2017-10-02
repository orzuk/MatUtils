% Compute distribution of selection coefficients observed in the population
% by weighting alleles according to their absorption time
% 
% Input: 
% s_birth_bins - values of selection coefficient for newly born mutations
% s_birth_hist - distribution of selection coefficient for newly born mutations
% N  - effective population size
%
% Output: 
% s_bins - values of observed selection distribution in the population
% s_hist - distribution of observed selection distribution in the population 
%
function [s_bins s_hist] = observed_selection_distribution(s_birth_bins, s_birth_hist, N)

theta = 0.001; % set some default mutation rate (doesn't affect results) 

num_s = length(s_birth_hist); s_hist = zeros(num_s, 1); 

s_bins = s_birth_bins;
for i=1:num_s % loop on all different s
    s_hist(i) = s_birth_hist(i) * absorption_time_by_selection(s_bins(i), theta, N);
end    

s_hist = normalize_hist(s_bins, vec2row(s_hist)); 




