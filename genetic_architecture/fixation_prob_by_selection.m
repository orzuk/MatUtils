% Compute mean time until absorbtion for allele with selection coefficient s
% The function computes the following expression:
% \int_{x_min}^{x_max} x^{weight_flag} * t_s(x) dx
% There is NO Normalization
%
% Input:
% s - selection coefficient (should be NEGATIVE for deleterious alleles)
% N - effective population size
% f_init - initial allele frequency (default is ~0, or 1/2N)
% x_max - maximal allele frequency (default is ~1, or 1-1/2N)
% weight_flag - how do we weight different allele frequencies.
%
% Output:
% P_fix - fixation probability by selection
%
function P_fix = fixation_prob_by_selection(s, N, f_init)

S = 4.*N.*s; % small s to big S
if(~exist('f_init', 'var') || isempty(f_init))
    f_init = 1/(2*N);
end

%P_fix = (1-exp(-2*s)) / (1-exp(-S)); % formula for f_init = 1/2N
P_fix = (1-exp(-S*f_init)) / (1-exp(-S)); % formula for f_init = 1/2N

