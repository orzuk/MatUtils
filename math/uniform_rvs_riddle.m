% Simulations used for riddle on uniform random variables
a = 5; b = 10; 
iters = 100000; 

% First part: generate uniform [0,1] r.v.s. from 3 uniform [a,b] r.v.s. where a and b are unknown
r = random('unif', a, b, 3, iters); % generate uniform r.v.s. in [a,b]
m = min(r); M = max(r); med = median(r);
u = (med - m) ./ (M - m); % transfer u to [0,1];
figure; hist(u, 100); 

% Second part: for which n can we have n (dependent!) uniform r.v.s. such
% that their sum is constant. For which k can we have each k of them be independent? 
% Note: The sum must be n/2
% Note: we can divide them into pairs, where each pair is x and 1-x and the
% pairs are independent. So any even n works



% Third part: Z = f(Y_1,..,Y_n) where the Y's are the order-statistics 
% and the Z's are uniformly distributed (though not independent) and Z is
% deterministic. Z can only be a permutation (??) 
n = 3;
r = random('unif', a, b, n, iters); % generate uniform r.v.s. in [a,b]

r = r .* (n/2) ./ repmat(sum(r), n, 1);
figure; hist(r(:), 100); 



