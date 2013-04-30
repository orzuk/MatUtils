% Compute log-likelihood of data vector according to a given MoG model. 
%
% Input Parameters:
% x - input data matrix 
% Mu - the means vector
% Sigma - the covariance matrix
% P - the probs. of each gaussian

% Output Parameters:
% LL - the log-likelihood of data 
%
% Written by Or Zuk 5/2007
%
function LL = MixtureOfGaussiansGetLikelihood(x, P, Mu, Sigma) 

num_gaussians = length(Sigma);
N=length(x); % # data points

p = zeros(N, num_gaussians);
for m=1:num_gaussians   
   p(:,m)= -((x'-Mu(m)).^2 ./ (2*Sigma(m)^2)) - log(Sigma(m)) - 0.5 * log(2*pi);   
end

LL = sum(exp(p) .* repmat(P, N, 1), 2);

LL = sum(log(LL));

