% Compute the mean variance explained by a locus with effect size beta
%
% Input:
% beta - effect size
% s_bins - vector of selection coefficients
% beta_bins - vector of effect sizes
% g_hist - joint dsitribution of beta and s
% N - effective population size
% x_vec - vector of different allele frequencies
%
% Output:
% V - total mean variance explained by one locus
% V_x_vec - density of mean variance explained according to derived allele frequency
%
function [V V_x_vec] = effect_size_to_mean_variance_explained(beta, s_bins, beta_bins, g_hist, N, x_vec)


num_beta = length(beta); % allow input vector
V = zeros(num_beta,1);

if(exist('x_vec', 'var') && (~isempty(x_vec))) % compute also variance explained by allele frequency
    V_x_vec = zeros(num_beta,length(x_vec));
end
for i_beta = 1:num_beta
    if(mod(i_beta,10)==0)
        run_beta = i_beta
    end
    [~, beta_ind] = min(abs(beta_bins-beta(i_beta))); % find the correct beta
    g_s_given_beta_hist = normalize_hist(s_bins, vec2row(g_hist(:,beta_ind))); % get conditional distribution of s for a given beta
    
    S = 4.*N.*s_bins; % convert to big S
    V(i_beta) = beta(i_beta)^2  * integral_hist(s_bins, g_s_given_beta_hist .* ...
        min(0.5, 1 ./ (1-exp(S)) + 1 ./ S) );
    
    if(exist('x_vec', 'var') && (~isempty(x_vec))) % compute also variance explained by allele frequency
        two_side_flag = 0; scale_mode = 'linear'; 
        for i=1:length(x_vec)
            V_x_vec(i_beta,i) = beta(i_beta)^2 * integral_hist(s_bins, g_s_given_beta_hist .* ...
                allele_freq_spectrum(x_vec(i), s_bins, N, two_side_flag, scale_mode, 1) );
        end
    end
end % loop on different betas
