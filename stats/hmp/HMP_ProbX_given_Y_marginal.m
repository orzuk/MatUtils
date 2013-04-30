% Compute P_X_i|Y, the prob. of hidden state X given observation Y 
% This is done for each of the 2^n combinations of X, and for a given 
% N, epsilon and p (one of the latter can be  given as a vector)
% The method is like the maple one (which was symbollicaly). Enumerate all
% 2^N configurations and sum the contribution P(X,Y).
% We should improve it to be less involved computationally (for each Y use
% dynamic programming/forward algorithm to compute P(Y))
function P_X_given_Y_marginal = HMP_ProbX_given_Y_marginal(eps_vec, p_vec, N, Y)


% One choice: Exhaustive search
P_X_given_Y = HMP_ProbX_given_Y(eps_vec, p_vec, N, Y);

p_len = length(p_vec); eps_len = length(eps_vec); 
P_X_given_Y_marginal = zeros(p_len, 2*N);

for i=1:2^N
    for j=1:N
        b = bitget(i-1, j)+1;
        P_X_given_Y_marginal(:,(j-1)*2+b) =  P_X_given_Y_marginal(:,(j-1)*2+b) + P_X_given_Y(:,i);
    end
end



% Other choice: Run Forward-Backward (much faster)


% % % % We assume initial condition : X_1 = 1
% % % p_len = length(p_vec);
% % % 
% % % P_X_given_Y = zeros(p_len, 2^N); % Make an array of Y probabilities
% % % 
% % % for X = 0:2^N-1  % Go over all the X's 
% % %     neq_bonds = sum(bitget(bitxor(X, bitshift(X,-1)), 1:N-1));
% % %     P_X = 0.5*((1-p_vec) .^ (N-1-neq_bonds) .* p_vec .^ neq_bonds)'; 
% % %         
% % %     neq_noise_bonds = sum(bitget(bitxor(X,Y), 1:N));
% % %     P_X_given_Y(:,X+1) = (1-eps_vec) .^ (N-neq_noise_bonds) .* eps_vec .^ neq_noise_bonds .* P_X';
% % %    
% % % end
% % % 
% % % P_X_given_Y = P_X_given_Y ./ repmat(sum(P_X_given_Y,2), 1, 2^N);


    
    
    


