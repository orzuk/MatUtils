
function [P_Ridge, Q_Ridge] = KroneckerRidge(Z_mat, options)

lambda = options.lambda;
[U, D, V] = svd(Z_mat);
[q, p] = size(Z_mat);
r = rank(Z_mat); %pp = p;
eig_vals = diag(D); % closed-form solution from Tibshirani's paper

% Compute beta, theta with penalties
c1 = -4*lambda*p^2;
c2 = 32*lambda^2*p + eig_vals.^4 .* (q-p);
%c3 = 4*lambda* (eig_vals.^4 - 16*lambda^2);

% Compute for increased numerical stability:
Discriminant = 64*lambda^2*p^2*eig_vals.^4 + (q-p)^2.*eig_vals.^8 + 64*lambda^2*p*(q-p)*eig_vals.^4;
% sqr1 = c2(1:r).^2 - 4.*c1.*c3(1:r)


beta = repmat(2*sqrt(lambda/p), q, 1);
theta = repmat(2*sqrt(lambda/q), p, 1);

if (p == q)
    beta(1:r) = sqrt( (eig_vals.^2 + 4*lambda)./p );
    theta=beta;
else
    beta(1:r) = real(sqrt( (-c2(1:r) - sqrt(Discriminant))  ./ (2.*c1) ));
    theta(1:r) = eig_vals(1:r).^2 .* beta(1:r) ./ (p .* beta(1:r).^2 - 4*lambda);
    theta(p .* beta(1:r).^2 - 4*lambda == 0) = 2*sqrt(lambda/q); % fix eigenvalues near zero
    theta(1:r) = max(theta(1:r), 2*sqrt(lambda/q));
end
P_Ridge = V*diag(theta)*V';
Q_Ridge = U*diag(beta)*U';
