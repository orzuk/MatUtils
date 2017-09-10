% Compute conditional probability of selection coefficient given frequency x in the population 
% 
% Input: 
% s - selection coefficient (can be a vector) 
% x_vec - allele frequency  
% s_bins - values of selection coefficient
% s_hist - density of selection coefficient
% N - effective population size 
% 
% Output: 
% p_s_given_x - Probability that selection coefficient is s for a random polymorphic allele with frequency x 
% p_x - Probability of seeing allele frequency x
% 
function [p_s_given_x p_x]= prob_selection_given_allele_freq(s, x_vec, s_bins, s_hist, N)

num_s = length(s); num_x = length(x_vec); 

p_s_given_x = zeros(num_s, num_x); % matrix of conditional probabilities
%p_x = zeros(num_s, num_x); % matrix of joint probabilities 

two_side_flag = 0; % look at derived allele frequency 

p_s_denominator = zeros(1, num_x); 
x_hist = zeros(length(s_bins), num_x);
for i=1:length(s_bins)
    x_hist(i,:) = x_vec .* exp(allele_freq_spectrum(x_vec, -s_bins(i), N, two_side_flag, 'log'));
    p_s_denominator = p_s_denominator + ...
        s_hist(i) .*  x_hist(i,:); % x_vec .* exp(allele_freq_spectrum(x_vec, s_hist(i), N, two_side_flag, 'log'));
%    p_x = p_x + s_hist(i) .* x_hist(i,:);
end
p_x = normalize_hist(x_vec, p_s_denominator);
p_x = repmat(p_x, num_s, 1);
% for j=1:num_x
%     p_s_denominator(i) = integral_hist(s_bins, s_hist .* x_hist(:,j)'); 
% end

for i=1:num_s % loop on selection coefficients
    [~, s_ind] = min(abs(s(i) - s_bins));
%    p_x(i,:) = s_hist(s_ind) .* x_vec .* exp(allele_freq_spectrum(x_vec, s(i), N, two_side_flag, 'log')); % multiply by allele freq !!!
    p_s_given_x(i,:) =  s_hist(s_ind) .* x_hist(s_ind,:) ./ p_s_denominator;
end

