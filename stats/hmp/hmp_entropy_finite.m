% Compute finite HMP entropy H_N for a given N, epsilon and p. 
% (The latter given as a vector)
% The method is like the maple one (which was symbollicaly). Enumerate all
% 2^N configurations and sum the p * log p expansion
function H_N = hmp_entropy_finite(eps, p_vec, N)


% We assume initial condition : X_1 = 1
p_len = length(p_vec);
H_N = zeros(p_len, 1);

P_Y = zeros(p_len, 2^N); % Make an array of Y probabilities

for X = 0:2^N-1  % Go over all the X's 
    neq_bonds = sum(bitget(bitxor(X, bitshift(X,-1)), 1:N-1));
    P_X = 0.5*((1-p_vec) .^ (N-1-neq_bonds) .* p_vec .^ neq_bonds)'; 
        
    for Y = 0:2^(N-1)-1  % Go over all the Y's. We can save half here due to symmetry !!!!!
        neq_noise_bonds = sum(bitget(bitxor(X,Y), 1:N));
        P_Y(:,Y+1) = (1-eps) .^ (N-neq_noise_bonds) .* eps .^ neq_noise_bonds .* P_X + P_Y(:,Y+1);
    end
   
end

% Now compute the entropy 
for Y = 0:2^(N-1)-1
    H_N = H_N - P_Y(:,Y+1) .* log(P_Y(:,Y+1));
end

% Transfer to binary logarithm
H_N = 2 .* H_N ./ log(2.0);
    

    
    
    


