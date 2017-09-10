% Compute probability of observing k out of n for a polymorphic allele
% 
% Input: 
% k - number of observed derived allele carriers 
% n - sample size
% s_bins - values of selection coefficient
% s_hist - density of selection coefficient
% N - effective population size 
% 
% Output: 
% p_k_n - Probability that a random polymorphic allele will have k derived allele carriers out of n
%
function p_k_n = prob_observed_allele_freq(k, n, s_bins, s_hist, N)

two_side_flag = 0; 
num_s = length(s_bins); num_k = length(k);
x_vec = (1:(2*N-1)) ./ (2*N); % set vector of derived allele frequencies 
p_k_n = zeros(num_k,1); cur_p_k_n = zeros(num_k,1); p_k_n_denominator = 0; 
for i=1:num_s % integrate over s 
    x_hist = exp(allele_freq_spectrum(x_vec, -s_bins(i), N, two_side_flag, 'log'));  % Compute allele freq distribution at equilibrium 
    for j=1:num_k
        cur_p_k_n(j) = ... % p_k_n(j) + s_hist(i) * ...
            integral_hist(x_vec, x_hist .* x_vec .^ k(j) .* (1-x_vec) .^ (n-k(j)));
    end
    p_k_n_denominator = 0 + ... % p_k_n_denominator + ... % s_hist(i) * ...
        integral_hist(x_vec, x_hist .* (1 - x_vec.^n  - (1-x_vec).^n));
    p_k_n = p_k_n + cur_p_k_n .* s_hist(i) ./ p_k_n_denominator; 
end
p_k_n = exp( vec2column(nchoosekln(n, k)) + log(p_k_n)); %  - log(p_k_n_denominator) ); % divide 
p_k_n = p_k_n ./ sum(s_hist); % we assume that the s bins are evenely spaced 
