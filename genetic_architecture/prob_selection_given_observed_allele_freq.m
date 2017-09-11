% Compute probability of selection coefficient s given that we observe k out of n for a polymorphic allele
% 
% Input: 
% s - selection coefficient (can be a vector) 
% k - number of observed derived allele carriers 
% n - sample size
% s_bins - values of selection coefficient
% s_hist - density of selection coefficient
% N - effective population size 
% 
% Output: 
% p_s_given_k_n - Probability that selection coefficient is s for a random polymorphic allele with k derived allele carriers out of n
%
function p_s_given_k_n = prob_selection_given_observed_allele_freq(s, k, n, s_bins, s_hist, N)

num_s = length(s); p_s_given_k_n = zeros(num_s, 1); 
two_side_flag = 0;
for i=1:num_s % allow a vector s as input 
    [~, cur_s_ind] = min(abs(s_bins - s(i)));
    p_s_given_k_n(i) = s_hist(cur_s_ind) * prob_observed_allele_freq(k, n, s(i), 1, N); % compute probability of this s and k and n 
    % Do we miss here the density of s? 
%     allele_freq_spectrum(x_vec, s, N, two_side_flag)
end

p_s_given_k_n = p_s_given_k_n ./ prob_observed_allele_freq(k, n, s_bins, s_hist, N); % normalize by phi(k, n) 
