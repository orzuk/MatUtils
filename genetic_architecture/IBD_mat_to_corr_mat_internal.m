% Convert binary IBD values to Gaussian correlations
%
% Input:
% IBD_mat - matrix of IBD sharing (probabilities between zero and one)
% f_vec - vector of allele frequencies
%
% Output:
% corr_mat -
% x_f_vec -
%
function [corr_mat x_f_vec] = IBD_mat_to_corr_mat_internal(IBD_mat, f_vec)

n = length(IBD_mat);

%prob11_mat = IBD_mat .*

corr_mat = zeros(n);
x_f_vec = norminv(f_vec); % get thresholds
for k=1:1 % loop on frequencies
    for i=1:n
        for j=i+1:n
            corr_mat(i,j) = fminbnd(@(x) mvncdf([x x]) - f_vec(k) * ...
                (IBD_mat(i,j) + f_vec(k)*(1-IBD_mat(i,j))), 0, 1);
        end
    end
end

