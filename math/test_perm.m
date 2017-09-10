% Test if permutation obeys certain conditions
% works only for small n (up to ~ 10) 
% Sampling can work for larger n but estimated prob. is very low 
function [good_frac cov_mat] = test_perm(n)

P = perms(1:n); % generate all possible permutations

good_vecs = double(mod( P(:,2:end) - P(:,1:end-1), n) - n/2 < 0);  

cov_mat = cov(good_vecs); 

good_vec = ones(factorial(n),1); 
for i=1:n-1 % Test 
    good_vec = good_vec .* double(mod( P(:,i+1) - P(:,i), n) - n/2 < 0);  
end

good_frac = sum(good_vec) / factorial(n); 

P_s_2_given_s_1_s_3 = sum(prod(good_vecs(:,[1:3]),2)) / sum(prod(good_vecs(:,[1 3]),2));
P_s_2_given_s_1_s_3_s4 = sum(prod(good_vecs(:,[1:4]),2)) / sum(prod(good_vecs(:,[1 3 4]),2));
