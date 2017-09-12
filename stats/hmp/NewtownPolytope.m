% Compute  Newton Polytope of inf. func. given an observation Y.
% For the binary-symmetric HMP (works only for two parameters)
function NP = NewtownPolytope(N,Y) % compute the Newton Polytope

S = zeros(N+1); % a sparse matrix

for X=0:2^N-1  % loop on all X's
    neq_bonds = sum(bitget(bitxor(X, bitshift(X,-1)), 1:N-1));
    %    P_X = 0.5*((1-p_vec) .^ (N-1-neq_bonds) .* p_vec .^ neq_bonds)';
    neq_noise_bonds = sum(bitget(bitxor(X,Y), 1:N));
    %   P_X_given_Y(:,X+1) = (1-eps_vec) .^ (N-neq_noise_bonds) .* eps_vec .^ neq_noise_bonds .* P_X';
    S(neq_bonds+1, neq_noise_bonds+1) = 1;

end

[x,y] = find(S);
I = convhull(x-1,y-1);
NP = [x(I) y(I)];

