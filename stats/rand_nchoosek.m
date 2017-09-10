% Pick at random k indices out of n
% 
% Input: 
% n - size of set to choose from
% k - size of set chosen
% 
% Output: 
% ind_vec - indices of vectors chosen
% 
function ind_vec = rand_nchoosek(n, k)

p = randperm(n);
ind_vec = zeros(n,1);
ind_vec(p(1:k)) = 1;

