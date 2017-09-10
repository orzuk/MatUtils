% Build the genetic relationship matrix according to Visscher's paper
% 
% Input: 
% snp_mat - matrix of SNPs for a set of individuals
% snp_type - binary (default) or diploid
% f_vec - allele frequencies 
% 
% Output: 
% G - matrix representing genetic relationship between pairs of individuals
% 
function G = genetic_relationship_matrix(snp_mat, snp_type, f_vec) %Build Genetic relashionship matrix

if(~exist('snp_type', 'var') || isempty(snp_type))
    snp_type = 'binary';
end

[n, m] = size(snp_mat); % m - # individuals, n - # snps
G = zeros(m);
f_mat = repmat(f_vec, 1, m);
switch snp_type
    case 'binary'
        G = (1/n) .* (snp_mat- f_mat)' * ((snp_mat- f_mat) ./ ...
            (f_mat .* (1-f_mat))); % Compute all non-diagonal elements
        % Need to deal with diagonal elements separately, but it doesn't
        % affect the regression
    case 'diploid'
        G = (1/n) .* (snp_mat- 2.*f_mat) * (snp_mat- 2.*f_mat)' ./ ...
            (2 .* f_mat .* (1-f_mat)); % Compute all non-diagonal elements
end




