% Compute the probability of Y according to an HMM
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary.
% 
% Input: 
% p - flip probability of markov chain
% eps - flip probability of channel
% Y - vector of observations
% 
% Output: 
% prob - probability of Y
% 
function prob = compute_HMM_condprob(p, eps, Y)

n = length(Y); % the vector length


% find the stationary distribution 
[V, D] = eig(p');

[val ind] = max(diag(D));

mu = V(:,ind) ./ sum(V(:,ind));

% Init the probabilities
p_y_when_x_1 = mu(1) * eps(1,Y(1));
p_y_when_x_2 = mu(2) * eps(2,Y(1));


normalize = 0; min_eps = 0.0000000001;

% Do the forward algorithm
for i=2:n-1
    new_p_y_when_x_1 = (p_y_when_x_1 * p(1,1) + p_y_when_x_2 * p(2,1)) * eps(1,Y(i));
    new_p_y_when_x_2 = (p_y_when_x_1 * p(1,2) + p_y_when_x_2 * p(2,2)) * eps(2,Y(i));
    
    p_y_when_x_1 = new_p_y_when_x_1;
    p_y_when_x_2 = new_p_y_when_x_2;
    
    % normalize ..
    if(min(p_y_when_x_1, p_y_when_x_2) < min_eps)
        p_y_when_x_1 = p_y_when_x_1 / min_eps;
        p_y_when_x_2 = p_y_when_x_2 / min_eps;
        normalize = normalize + log(min_eps);
    end
        
end


% Now do one final iteration. We must not have normalization here !!!! 
final_p_y_when_x_1 = (p_y_when_x_1 * p(1,1) + p_y_when_x_2 * p(2,1)) * eps(1,Y(i));
final_p_y_when_x_2 = (p_y_when_x_1 * p(1,2) + p_y_when_x_2 * p(2,2)) * eps(2,Y(i));

prob = (final_p_y_when_x_1 + final_p_y_when_x_2)/(p_y_when_x_1 + p_y_when_x_2);
