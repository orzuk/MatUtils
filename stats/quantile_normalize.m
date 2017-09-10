% Perform quantile normalization on a matrix.
% Sort each column and replace each value by the mean 
% of all values with the same column ranks
function D_norm = quantile_normalize(D)

[m n] = size(D); 
[D_norm D_inds] = sort(D); % sort columns 

vals_vec = mean(D_norm,2); 
for i=1:n
    D_inds(:,i) = inv_perm(D_inds(:,i));
end
D_norm = vals_vec(D_inds); % put the right inds 


