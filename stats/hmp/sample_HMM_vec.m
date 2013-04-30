% Sample a vector of observation from a given HMM 
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary
function Y = sample_HMM_vec(p, eps, sample_len, num_samples)


% find the stationary distribution 
[V, D] = eig(p');
[val ind] = max(diag(D));
mu = V(:,ind) ./ sum(V(:,ind));


X = ones(num_samples, sample_len); Y = ones(num_samples, sample_len);

% randomize x
rand_vec = rand(num_samples, sample_len);

X(:,1) = (rand_vec(:,1) < mu(1)) + 1;
for i=2:sample_len
    two_ind = find(rand_vec(:,i) > p(X(:,i-1), 1));
    X(two_ind,i) = 2;
    
% % % %      if(rand_vec(i) > p(X(i-1), 1))
% % % %          X(i) = 2;
% % % %      end
     
% % %       X(i) = 1 + ceil(rand_vec(i) > p(X(i-1), 1));

end

% randomize y
rand_vec = rand(num_samples, sample_len);

for i=1:sample_len
    two_ind = find(rand_vec(:,i) > eps(X(:,i), 1));
    Y(two_ind,i) =2;
end

