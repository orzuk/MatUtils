% Compute prob. of i.i.d. state to be equal given dirty observation
% p is the ones probability for x.
% eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary.
% The function returns the phi value, which is : 
% phi = (1/N) log Pr(All M chains are equal in the first N bases)
%
function equal_log_prob = compute_iid_eq_prob_given_eq_obser(p, eps)



%%%equal_log_prob = (2*p.^2 -2*p + 1)*(2*eps.^2-2*eps+1) ./ (  2*eps.*(eps-1) .* (2*p-1).^2 +  2*p.^2 -2*p + 1);


equal_log_prob = (2*p.^2 -2*p + 1)*(2*eps.^2-2*eps+1) ./ ( (2*p.^2 -2*p + 1)*(2*eps.^2-2*eps+1) + 4*p.*(1-p).*eps.*(1-eps)  );

%%equal_log_prob = log2(equal_log_prob);


