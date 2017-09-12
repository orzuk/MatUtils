% Computes P_Y, the probability of the observation Y. 
% This is done for each of the 2^n combinations of Y, and for a given 
% N, epsilon and p (one of the latter can be  given as a vector)
% The method is like the maple one (which was symbollicaly). Enumerate all
% 2^N configurations and sum the contribution P(X,Y).
% We should improve it to be less involved computationally (for each Y use
% dynamic programming/forward algorithm to compute P(Y))
%
function P_Y = HMP_ProbY_finite(eps_vec, p_vec, N)


% We assume initial condition : X_1 = 1
p_len = length(p_vec);

P_Y = zeros(p_len, 2^N); % Make an array of Y probabilities

for X = 0:2^N-1  % Go over all the X's 
    neq_bonds = sum(bitget(bitxor(X, bitshift(X,-1)), 1:N-1));
    P_X = 0.5*((1-p_vec) .^ (N-1-neq_bonds) .* p_vec .^ neq_bonds)'; 
        
    for Y = 0:2^N-1  % Go over all the Y's. We can save half here due to symmetry !!!!!
        neq_noise_bonds = sum(bitget(bitxor(X,Y), 1:N));
        P_Y(:,Y+1) = ((1-eps_vec) .^ (N-neq_noise_bonds) .* eps_vec .^ neq_noise_bonds .* P_X' + P_Y(:,Y+1)')';
    end
   
end


    
    
    


