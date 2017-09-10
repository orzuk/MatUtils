% Compute joint prob. dist. for family members
%
% Input:
% f_vec - marginal frequencies of SNPs
% kinship - kinship coefficient
% U - probability of sharing 0,1,2 alleles (needed for genotype mode)
% mode - use genpotypes or alleles (assumes two alleles are independent) 
% 
% Output:
% T - matrix of joint genotypes/allele frequencies. T(i,j) = Pr(x=i, x_R=j)
%
function T = joint_familial_marginal_probs(f_vec, kinship, U, mode, varargin)

if(~exist('mode', 'var') || isempty(mode)) % default: compute allele transistion probabilities
    mode = 'allele';
end
N = length(f_vec);
switch mode
    case 'allele'
        T = zeros(2,2,N);
        for i=0:1
            for j=0:1
                T(i+1,j+1,:) = f_vec.^i .* (1-f_vec).^(1-i)  .* ...
                    (kinship .* (1-(i-j)^2) + (1-kinship) .* f_vec.^j .* (1-f_vec).^(1-j));
            end
        end        
    case 'genotype'
        T = zeros(3,3,N);
        T(1,1,:) = 1-f_vec; T(1,2,:) = f_vec; % set manually transition probabilities when one allele is shared
        T(2,1,:) = 0.5.*(1-f_vec); T(2,2,:) = 0.5; T(2,3,:) = 0.5.*f_vec;
        T(3,2,:) = 1-f_vec; T(3,3,:) = f_vec; 
        for i=0:2 % AA, AB, BB
            for j=0:2
                T(i+1,j+1,:) = reshape(f_vec.^i .* (1-f_vec).^(2-i) .* 2^(i*(2-i)), 1, 1, N) .* ...
                    (U(3) .* (i == j) + U(2) .* T(i+1,j+1,:) + ...
                    U(1) .* reshape(f_vec.^j .* (1-f_vec).^(2-j) .* 2^(j*(2-j)), 1, 1, N));
            end
        end
end
