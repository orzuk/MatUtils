% Compute the intersection state (X) probability given observation (Y)
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary
function Y_equal_prob = comp_prob_X_equal_given_Y_equal(p, eps, HMM_flag)

% Help variables
p2 = 2*p.^2-2*p+1;
eps2 = 2*eps.^2-2*eps+1;

if(HMM_flag == 0) % Here the i.i.d. case

    Y_equal_prob = (p2 .* eps2) ./ (p2 + 2*eps.*(eps-1).*(2*p-1).^2);
else  % Here do the HMM
    
    % First the nominator
    Y_equal_prob = p2 .* (2*eps-1).^2 + sqrt( p2.^2 -8*eps.*(1-eps).*eps2.*(2*p-1).^2 ); 
    Y_equal_prob = Y_equal_prob ./ ( Y_equal_prob + 8*p.*eps.*(1-p).*(1-eps) );
    
end