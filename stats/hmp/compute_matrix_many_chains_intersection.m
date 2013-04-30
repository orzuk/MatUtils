% Compute prob. of intersection for different HMP realizations 
% p is the transition probability matrix S*S, where we have M different 
% markov chains and S = [M/2]+1
% Currently we use only odd chains for simplicity
%
% eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary.
% The function returns the phi value, which is : 
% phi = (1/N) log Pr(All M chains are equal in the first N bases)
%
function intersection_log_prob = compute_matrix_many_chains_intersection(p, eps, M)



% Prepare a binomial table 
binom_tab = zeros(M+1, M+1);


for i=1:M+1
    binom_tab(i,1:i) = binom(i-1, [0:i-1]);
end
binom_tab;
eps_power_vec = eps .^ [0:M];
one_minus_eps_power_vec = (1-eps) .^ [0:M];
p_power_vec = p.^ [0:M];
one_minus_p_power_vec = (1-p) .^ [0:M];


S = floor(M/2)+1;    % Set the matrix dimension

A = zeros(S,S);

% Compute the matrix elements one-by-one
% Note : State i means that the minority has exactly i-1 elements !!
for j=0:S-1
    Y_given_X = (eps*(1-eps))^j*(eps^(M-2*j)+(1-eps)^(M-2*j));
    for i=0:S-1
        % First i->j
        k_vec = [max(0,i-j):min(i,M-j)];
        X_given_prev_X = sum(binom_tab(i+1,k_vec+1) .* binom_tab(M-i+1,j-i+k_vec+1) .* p_power_vec(j-i+2*k_vec+1) .* one_minus_p_power_vec(M+i-j-2*k_vec+1));
        % Now i->M-j
        k_vec = [max(0,i-M+j):min(i,j)];
        X_given_prev_X = X_given_prev_X + sum(binom_tab(i+1,k_vec+1) .* binom_tab(M-i+1,M-j-i+k_vec+1) .* p_power_vec(M-j-i+2*k_vec+1) .* one_minus_p_power_vec(i+j-2*k_vec+1));
        A(i+1,j+1) = Y_given_X * X_given_prev_X;
    end
end

% Do boundary correction for S=M/2 if M is even
if(mod(M,2)==0)
    A(:,S) = 0.5*A(:,S);        
end

lambda = eig(A);


intersection_log_prob = log2(max(abs(lambda)));


