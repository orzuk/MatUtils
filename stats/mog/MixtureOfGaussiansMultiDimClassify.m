% Classfy a vector of points to the best fitting multi-variate Gaussians.  
% For each data point give a class according to which gaussian does it
% belong out of the multi-dimensional MOG. The gaussian that best explain
% the point is chosen, i.e. the g maximizing: Pr(g|x)
%
% The method is simply calculating the Mahalahonis distance.
%
% Input:
% x - input data matrix 
% P_vec - the probs. of each Gaussian
% Mu_vecs - the means vector
% Sigma_mats - the covariance matrix
%
% Output:
% C - the chosen classes covariance matrix
% CondProbs - the conditional probabilities for each data point to be in
% each class - Currently not used !!! 
%
% Written by Or Zuk 5/2007
%
function C = MixtureOfGaussiansMultiDimClassify(x, P_vec, Mu_vecs, Sigma_mats) 

num_gaussians = length(Sigma_mats);
dim=size(x,1); % Data dimension
N=size(x,2); % # data points

p = zeros(N, num_gaussians);
for m=1:num_gaussians
%%   p(:,m)=exp(- sum(((x'-repmat(Mu_vecs(m,:),N,1)) * inv(Sigma_mats{m})  .* ...
%%       (x'-repmat(Mu_vecs(m,:),N,1)   ))  ./ 2,2 ) ) ./ (sqrt((2*pi)^dim.*det(Sigma_mats{m})));
   p(:,m)=(- sum(((x'-repmat(Mu_vecs(m,:),N,1)) * inv(Sigma_mats{m})  .* ...
       (x'-repmat(Mu_vecs(m,:),N,1)   ))  ./ 2,2 ) ) - log(det(Sigma_mats{m})); % No need to keep the sqrt(2*pi) around
end

[DummyVal C] = max(p');
