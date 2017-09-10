% Compute covariance for two matrices 
function C = my_cov(X,Y)
n = size(X, 1)
C = (X' * Y) ./ (n-1) - mean(X)' * mean(Y) .* (n/(n-1));
