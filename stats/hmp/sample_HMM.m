% Sample an observation from a given HMM 
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary
function Y = sample_HMM(p, eps, sample_len)


% find the stationary distribution 
[V, D] = eig(p');
[val ind] = max(diag(D));
mu = V(:,ind) ./ sum(V(:,ind));


X = ones(1, sample_len); Y = ones(1, sample_len);

% randomize x
rand_vec = rand(1, sample_len);

X(1) = (rand_vec(1) < mu(1)) + 1;
for i=2:sample_len
     if(rand_vec(i) > p(X(i-1), 1))
         X(i) = 2;
     end
     
% % %       X(i) = 1 + ceil(rand_vec(i) > p(X(i-1), 1));

end

% randomize y
rand_vec = rand(1, sample_len);

for i=1:sample_len
    if(rand_vec(i) > eps(X(i), 1))
        Y(i) = 2;
    end
end

