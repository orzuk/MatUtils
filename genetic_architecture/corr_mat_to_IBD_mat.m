% Convert Gaussian correlations to binary IBD values 
%
% Input:
% corr_mat - Matrix of correlations of Gaussin random variables 
% f_vec - vector of allele frequencies
%
% Output:
% IBD_mat - Matrix of IBD sharing (probabilities between zero and one)
% x_f_vec - thresholds relating Gaussian and binary variables
%
function [IBD_mat x_f_vec] = corr_mat_to_IBD_mat(corr_mat, f_vec)

n = length(corr_mat);
if(length(f_vec) == 1) % same frequencies 
    f_vec = repmat(f_vec, 1); 
end
IBD_mat = zeros(n);
x_f_vec = norminv(f_vec); % get thresholds
for k=1:1 % loop on frequencies
    for i=1:n
        for j=i+1:n
            IBD_mat(i,j) = mvncdf([x_f_vec(i) x_f_vec(j)], [0 0], [1 corr_mat(i,j); corr_mat(i,j) 1]);
%            IBD_mat(i,j) = (tmp_IBD - f_vec(k).^2) / (f_vec(k)-f_vec(k).^2);
        end
    end
end
IBD_mat = IBD_mat + IBD_mat'; 
IBD_mat = (IBD_mat - f_vec' * f_vec) ./ sqrt( (f_vec .* (1-f_vec))' *  (f_vec .* (1-f_vec)) );  % correct 
IBD_mat = IBD_mat - diag(diag(IBD_mat)) + eye(n); % get symmetric and self correlations 


