% Compute fixation probability for newly born allele 0ith selection coefficient s
% The function computes the following expression:
% (1-e^(-2s)) / (1-e^(-4Ns))
%
% Input:
% s - selection coefficient (should be NEGATIVE for deleterious alleles)
% N - effective population size
% f_init - initial allele frequency (default is ~0, or 1/2N)
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
if(s == 0) % take limit (to avoid Nans)
    P_fix = 1./(2.*N);
else
    P_fix = (1-exp(-S*f_init)) / (1-exp(-S)); % formula for general f_init = 1/2N
end

