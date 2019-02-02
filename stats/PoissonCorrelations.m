% check pairwise correlations 
p = 0.01; 
epsilon = 0.0005; 
M = 10; % number of variables in each group

Sigma = eye(2*M)*p*(1-p);
Sigma(1:M, (M+1):(2*M)) = epsilon; 
Sigma((M+1):(2*M), 1:M) = epsilon; 

rho = epsilon / (p*(1-p))


isposdef(Sigma)
figure; imagesc(Sigma); colorbar;

n = 100000; 
X = mvnrnd(zeros(2*M, 1), Sigma, n);


S = corr(X)
figure; imagesc(S); colorbar;
median(S(1:M, (M+1):(2*M)))




corr(X(:,1), X(:, 2*M))

corr(sum(X(:,1:M), 2), sum(X(:, (M+1):(2*M)), 2))


