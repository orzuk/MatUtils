% Sample points from subspaces
% Input:
% p - dimension of entire space
% d - dimension of subspaces
% K - number of subspaces
% n - number of points
% sigma - noise level
% theta - angle between subspaces (optional)
%
% Output:
% X - data matrix
% B - basis matrix
%
function [X, B, theta_arr] = subspace_sample(p, d, K, n, sigma, theta)

% first sample basis
B = randn(p, K, d);
for i=1:K
    B(:,i,:) = reshape(orth(reshape(B(:,i,:), p, d)), p, 1, d);  % make orthonormal basis
end
theta_arr = zeros(K);
for j=1:(K-1) % compute their angles
    for k=(j+1):K
        theta_arr(j,k) = subspace(reshape(B(:,j,:), p, d), reshape(B(:,k,:), p, d));
    end
end
theta_mean = sum(theta_arr(:)) / nchoosek(K,2);

if(exist('theta', 'var') && (~isempty(theta))) % force angle between subspaces to be theta
    B0 = orth(randn(p, d)); % take standard directions (w.l.o.g)

    % Determine alpha:
    alpha = theta / theta_mean; % crude approximation
    for i=1:K
        B(:,i,:) = reshape(orth(alpha*reshape(B(:,i,:), p, d) + (1-alpha)*B0), p, 1, d);  % make orthonormal basis
    end
    for j=1:(K-1) % compute their angles again
        for k=(j+1):K
        theta_arr(j,k) = subspace(reshape(B(:,j,:), p, d), reshape(B(:,k,:), p, d));
        end
    end
end

% Next sample points
Z = randn(d, n*K); X = zeros(p, n);
for i=1:K
    X(:, ((i-1)*n/K+1):(i*n/K)) = reshape(B(:,i,:), p, d) * Z(:, ((i-1)*n/K+1):(i*n/K));
end
X = X + randn(p, n) .* sigma; % add noise








