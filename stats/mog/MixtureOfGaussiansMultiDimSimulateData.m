% Simulate samples from a Multi-Dimensional MOG distribution 
%
% Input:
% P_vec - probs. of each gaussian
% Mu_vecs - means vector
% Sigma_mats - covariance matrix
% n - number of data points to simulate 
%
% Output:
% x - data vec 
%
% Written by Or Zuk 5/2007
%
function x = MixtureOfGaussiansMultiDimSimulateData(P_vec, Mu_vecs, Sigma_mats, n) 

num_gaussians = length(Sigma_mats);
dim = length(Sigma_mats{1});

% First determine which gaussian each one gets
cum_P_vec = cumsum(P_vec);
G = rand(1,n); R = zeros(1,n);

for i=num_gaussians:-1:1
    R(find(G <= cum_P_vec(i))) = i; 
end

% Simulate the Gaussian data itself
G = randn(dim,n); x = zeros(dim,n);
for i=1:num_gaussians
    A  = sqrtm(Sigma_mats{i});
    x(:,find(R == i)) = A*G(:,find(R == i))  + repmat(Mu_vecs(i,:)', 1, sum(R == i));
end
