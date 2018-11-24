% Compute mean time until absorption for allele with selection coefficient
% s exactly using direct Markov chain approach
%
% Input:
% s - selection coefficient (should be POSITIVE for deleterious alleles)
% N - effective population size
%
% Output:
% T - mean time until absorption for a newly born allele (in # of generation) 
%
function T = absorption_time_by_selection_numeric(s, N)

% TEMP! get exact absorption probability to determine ratio: 
% TEMP: compute for small example also absorption probability at stationary
M = FisherWright_ComputeMarkovMatrix(N, s, 'exact', 1); % compute Markov matrix
mu=1; M(1,2)=mu; M(1,1)=1-mu; M(end,2)=mu; M(end,end)=1-mu; % add small mutation to make process ergodic

%[T_vec, T_mat] = MarkovChainAbsoptionTime(M, [1 2*N+1]); 
pi_vec = markov_chain_stationary_dist(M); % pi_poly_vec = [0 pi_vec(2:end-1)' ./ sum(pi_vec(2:end-1)) 0]'; 
T = pi_vec(1)+pi_vec(end); T = (1-T)/T; % T = pi_poly_vec' * M; T = T(1)+T(end)





