% Convert binary IBD values to Gaussian correlations
%
% Input:
% IBD_mat - Matrix of IBD sharing (probabilities between zero and one)
% f_vec - vector of allele frequencies
%
% Output:
% corr_mat - Matrix of correlations of Gaussin random variables
% x_f_vec - thresholds relating Gaussian and binary variables
%
function [corr_mat x_f_vec] = IBD_mat_to_corr_mat(IBD_mat, f_vec)

n = length(IBD_mat);
if(length(f_vec) == 1) % same frequencies
    f_vec = repmat(f_vec, 1);
end

%prob11_mat = IBD_mat .*
corr_mat = zeros(n);
x_f_vec = norminv(f_vec); % get thresholds
%for k=1:1 % loop on frequencies
for i=1:n
    for j=i+1:n
        corr_mat(i,j) = fminbnd(@(x) ...
            (mvncdf([x_f_vec(i) x_f_vec(j)], [0 0], [1 x; x 1]) - ...
            ( f_vec(i)*f_vec(j) + ...
            IBD_mat(i,j) * sqrt( f_vec(i)*f_vec(j)*(1-f_vec(i))*(1-f_vec(j))) )).^2, 0, 1);
        %             f_vec(k) * ...
        %                 (IBD_mat(i,j) + f_vec(k)*(1-IBD_mat(i,j)))).^2, 0, 1);
    end
end
%end
corr_mat = corr_mat + corr_mat' + eye(n); % get symmetric and self correlations





