% Compute effect sizes and selection coefficient density based on Eyre-Walker's
% distribution from Eyre-Walker, PNAS, 2010.
% http://www.pnas.org/content/107/suppl.1/1752.abstract
% The model is:
% s ~ Gamma(k, theta/k)
% beta = delta*s.^tau*(1+epsilon)
% epsilon ~ N(0, sigma_epsilon^2)
%
% Input:
% tau - power relating s and beta
% sigma_epsilon - noise in coupling of s and beta
% delta - constant relating s and beta (in the paper it is drawn uniformly as +/- 1 with prob. 1/2) 
% theta - scale parameter of the Gamma distribution
% k - shape parameter of the Gamma distribution
% N - effective population size 
% n - number of loci to sample
% link_function - coupling between s and beta (default 'linear', also possible 'sigmoid')
%
% Output:
% beta - list of sampled effect sizes
% s - list of sampled selection coefficients
%
function [beta s] = ...
    EyreWalker_dist_pdf(tau, sigma_epsilon, delta, ... % parameters relating s and beta
    theta, k, N, beta, s, link_function) % parametring specifying the gamma distribution for beta

if(~exist('N', 'var') || isempty(N))
    N = 10000; % effective population size of humans
end
if(~exist('link_function', 'var') || isempty(link_function)) % set default link function
    link_function = 'linear';
end

s = gamrnd(k, theta/k, n, 1); % simulate from a Gamma distribution
eps = randn(n, 1) .* sigma_epsilon; % randomize
%sign_vec = sign(rand(n, 1) - 0.5);

switch link_function
    case 'linear' % This is the original model in the paper 
        beta = abs(delta .* s.^tau .* (1+eps)); % .* sign_vec % compute affect size
    case 'sigmoid' % This is an alternative model 
        beta = (  logit(  min_abs(delta .* (s .* (1+eps) - 0*median(s)) ./ (max(s)), 0.99999999999999))); % effect size for
        beta = beta - min(beta);
end
s = s ./ (4*N); % move from big S to small s

