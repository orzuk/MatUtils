
function [P_Ridge, Q_Ridge] = KroneckerRidge(Z_mat, options)

lambda = options.lambda;
[U, D, V] = svd(Z_mat);
[q, p] = size(Z_mat);
r = rank(Z_mat); pp = p;
eig_vals = diag(D); % closed-form solution from Tibshirani's paper

% Compute beta, theta with penalties
c1 = -4*lambda*pp^2;
c2 = 32*lambda^2*pp + eig_vals.^4 .* (q-pp);
c3 = 4*lambda* (eig_vals.^4 - 16*lambda^2);

beta = repmat(2*sqrt(lambda/pp), q, 1);
beta(1:r) = real(sqrt( (-c2(1:r) - sqrt(c2(1:r).^2 - 4.*c1.*c3(1:r)))  ./ (2.*c1) ));
theta = repmat(2*sqrt(lambda/q), pp, 1);
theta(1:r) = eig_vals(1:r).^2 .* beta(1:r) ./ (pp .* beta(1:r).^2 - 4*lambda);
theta(pp .* beta(1:r).^2 - 4*lambda == 0) = 2*sqrt(lambda/q); % fix eigenvalues near zero
theta = max(theta(1:r), 2*sqrt(lambda/q));
P_Ridge = V*diag(theta)*V';
Q_Ridge = U*diag(beta)*U';
