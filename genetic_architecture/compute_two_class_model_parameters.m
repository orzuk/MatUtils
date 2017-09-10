% Compute all parameters of interest for two class model
%
% Input:
% s_null_vec - vector of possible selection coefficients for the null alleles
% f_rare_vec - thresholds for an allele to be considered rare
% alpha_vec - fraction of null alleles at birth
% rare_cumulative_per_gene - total expected number of rare alleles in gene
% N - effective population size
% output_file - where to save all computed variables
%
% Output:
% two_class_stat_struct - structure containing many additional parameters
% w_x_null_mat - Matrix of weighted time below freq. <= f for alleles as function of f and selection coefficient s.  Notation: int_{0}^{f(i)} x T_{s(j)}(x) dx.
% w_x_harmless - Vector of weighted time below freq. <= f for NEUTRAL alleles.  Notation: int_{0}^{f(i)} x T_0(x) dx.
% w_all - Vector of total time of allele being polymorphic. Notation: int_{0}^{1} x T_{s(j)}(x) dx
% xi_prob_leq_f - Matrix of probability of allele in mixture being below prob. f. Notation: \xi( <= f)
% frac_null_by_freq_cumulative - Prob(allele null | freq. <= f). Notation:  \alpha_{s^*}(<= f)
%
function [two_class_stat_struct ...
    w_x_null_mat w_x_harmless w_all xi_prob_leq_f frac_null_by_freq_cumulative] = ...
    compute_two_class_model_parameters(s_null_vec, f_rare_vec, alpha_vec, rare_cumulative_per_gene, N, output_file)

w_x_null_mat = zeros(length(s_null_vec), length(f_rare_vec)); % for each s and allele-freq cutoff f, compute absorption time
w_x_harmless = zeros(1, length(f_rare_vec)); % absorption time when selection coefficient is zero
for i=1:length(f_rare_vec) % loop on different rare-allele cutoffs
    if(mod(i, 10) == 0)
        compute_two_class_model_w_rare_allele_cutoff = i
    end
    w_x_harmless(i) = absorption_time_by_selection(0, 1, N, 0, f_rare_vec(i), 'freq');
    for j=1:length(s_null_vec) % loop on different selection coefficients
        w_x_null_mat(j,i) = absorption_time_by_selection(s_null_vec(j), 1, N, 0, f_rare_vec(i), 'freq');
    end
end
w_all = alpha_vec .* absorption_time_by_selection(s_null_vec, 1, N, 0, 0.9999999999999, 'freq') + ...
    (1-alpha_vec) .* absorption_time_by_selection(0, 1, N, 0, 0.99999999999999, 'freq');
xi_prob_leq_f = (alpha_vec) .* w_x_null_mat + (1-alpha_vec) .* repmat(w_x_harmless, length(s_null_vec), 1); % here count how many rare alleles we should have !!
xi_prob_leq_f = (xi_prob_leq_f ./ repmat(w_all, 1, length(f_rare_vec))) .* rare_cumulative_per_gene;
frac_null_by_freq_cumulative = alpha_vec .* w_x_null_mat ./ ...
    (alpha_vec .* w_x_null_mat + (1-alpha_vec) .* repmat(w_x_harmless, length(s_null_vec), 1));

% Compute auxillary parameters
% Compute mean and median for each s
num_s_null = length(s_null_vec);
mean_x = zeros(num_s_null,1);
median_x = zeros(num_s_null,1); median_mixture_x = zeros(num_s_null,1);
mean_het_x = zeros(num_s_null,1); total_het_x = zeros(num_s_null,1); 
integral_phi_x = zeros(num_s_null,1);
normalization_factor_x = zeros(num_s_null,1); 
S_vec = s_null_vec .* N .* 4;
for i=1:length(s_null_vec)
    S = s_null_vec(i) * N * 4;
    integral_phi_x(i) = phi_s_integral(0.999999999,-S, 0) - phi_s_integral(0.000000001,-S, 0);
    mean_x(i) = phi_s_integral(0.999999999,-S, 2) - phi_s_integral(0.000000001,-S, 2);        
    normalization_factor_x(i) = phi_s_integral(0.999999999,-S, 1) - phi_s_integral(0.000000001,-S, 1);    
    mean_x(i) = mean_x(i) ./ normalization_factor_x(i); % normalize to get mean of polymorphic alleles 
    median_x(i) = fzero(@(x) phi_s_integral(x,-S, 1) - ...
       0.5 .* phi_s_integral(0.999999999,-S, 1) - 0.5*phi_s_integral(0.000000001,-S, 1), mean_x(i)); % solution of phi(x)-phi(0) = (1/2)*( phi(1)-phi(0))
    mean_het_x(i) = phi_s_integral(0.999999999,-S, -2) - phi_s_integral(0.000000001,-S, -2);
    mean_het_x(i) = mean_het_x(i) ./  normalization_factor_x(i);
    total_het_x(i) = phi_s_integral(0.999999999,-S, 'var') - phi_s_integral(0.000000001,-S, 'var');
end
mean_mixture_x = frac_null_by_freq_cumulative(:,end) .* mean_x + ...
    (1-frac_null_by_freq_cumulative(:,end)) .* mean_x(1); % mean is just a weighted average

for i=1:length(s_null_vec)
    [~, median_mixture_x(i)] = min( (xi_prob_leq_f(i,:) - 0.5).^2); 
    median_mixture_x(i) = f_rare_vec(median_mixture_x(i));
%     S = s_null_vec(i) * N * 4;
%     median_mixture_x(i) = fzero(@(x) frac_null_by_freq_cumulative(i,end) .* (phi_s_integral(x,-S, 1) - ...
%         (0.5 .* phi_s_integral(0.999999999,-S, 1) - 0.5*phi_s_integral(0.000000001,-S, 1))) + ...
%         (1-frac_null_by_freq_cumulative(i,end)) .* (phi_s_integral(x, 0, 1) - ...
%         (0.5 .* phi_s_integral(0.999999999,0, 1) - 0.5*phi_s_integral(0.000000001,0, 1))), mean_x(i));
end

bayes_factor_vec = w_x_null_mat ./ ...
    repmat(w_x_harmless, length(s_null_vec), 1); %  enrichment of neutrals vs. nulls at different allele frequencies

two_class_stat_struct = struct('mean_x', mean_x, 'median_x', median_x, 'mean_het_x', mean_het_x, 'total_het_x', total_het_x, ...
    'integral_phi_x', integral_phi_x, 'normalization_factor_x', normalization_factor_x, ...
    'mean_mixture_x', mean_mixture_x, 'median_mixture_x', median_mixture_x, 'bayes_factor_vec', bayes_factor_vec);


% if(exist('output_file', 'var') && (~isempty(output_file))) % save output to file 
%     save(output_file, 'two_class_stat_struct', ...
%         'w_x_null_mat', 'w_x_harmless', 'w_all', 'xi_prob_leq_f', 'frac_null_by_freq_cumulative');
% end
