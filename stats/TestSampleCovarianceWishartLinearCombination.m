% Test linear combinations of elements of Covariance matrix which as a
% matrix has a Wishart distribution under Normality
p = 10; % size of two blocks
q = 10; 

Sigma11 =  posdefrnd(p); 
Sigma22 = posdefrnd(q); 

Sigma11 = diag(diag(Sigma11)); Sigma22 = diag(diag(Sigma22)); 

Sigma = [Sigma11, zeros(q, p); zeros(p, q), Sigma22]; d=p+q;
% Sigma = kron(Sigma11, Sigma22); d=p*q; 
  
iters = 10000; 

X = mvnrnd(zeros(d, 1), Sigma, iters); 

T_offdiag_sqr = zeros(iters, 1); % Compute off-diagonal 
for i=1:iters
   S = X(i,:)'*X(i,:);
   T_offdiag_sqr(i) = sum(sum( S((p+1):end, 1:q).^2 )); 
end    


figure; hist(T_offdiag_sqr, 100); 




