% An estimator for theta, the overall (population) mutation rate form a sample
% Input:
% X - a matrix of genotypes in the sample
% estimator_type - how to estimate theta.
%
% Output:
% Theta - an estimator for the overall mutation rate
%
function Theta = EstimateThetaCoalescent(X) % , estimator_type)

[n, L] = size(X); % get number of individuals and number of polymorphic sites

n_sites = length(find(max(X)-min(X))); % find polymorphic sites


n_alleles = size(unique(X, 'rows'), 1); % Count # alleles

%switch estimator_type
%    case 'MLE'
Theta(1) = fzero(@(t) internal_mean_allele_num(t, n) - n_alleles, n_alleles);
%    case 'het'
Theta(2) = sum(sum(X) .* (n-sum(X))) / (n*(n-1)/2);  % average pairwise differences
%    case 'num_alleles'
Theta(3) = n_sites / sum(1 ./ (1:(n-1))); % count total number of sites
%   end

% Compute expected number of distinct alleles in infinite alleles model
% Output:
% theta - mutation rate
% n - sample size 
%
% Output:
% a - expected number of alleles
%
function a = internal_mean_allele_num(theta, n)

a = 0;
for i=1:n % use scalar version to allow optimziation for vector theta
    a = a + theta ./ (theta + i-1);
end


