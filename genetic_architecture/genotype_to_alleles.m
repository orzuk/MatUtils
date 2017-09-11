% Expand genotype to pairs of alleles
% Example: ...
%
% Input:
% g_vec - vector of genotypes (0,1,2 in each)
% g_probs - probability of each genotype
%
% Output:
% x_vec - vector of alleles (in adjacent pairs)
% x_probs - probability of each allele
%
function [x_vec x_probs] = genotype_to_alleles(g_vec, g_probs)

[m N] = size(g_vec);


num_allele_vecs = prod(2.^(g_vec == 1),2);
x_vec = zeros(sum(num_allele_vecs), 2*N); x_probs = zeros(sum(num_allele_vecs),1);

tmp_all_x_vec =  zeros(max(num_allele_vecs), 2*log2(max(num_allele_vecs)));
tmp_all_x_vec(:,1:2:end) = dec2bin(0:max(num_allele_vecs)-1) - '0'; % take all 2^N possibilities
tmp_all_x_vec(:,2:2:end) = 1-tmp_all_x_vec(:,1:2:end);

ctr=1;
for i=1:m % loop over all genotypes (can't do it in matrix form for now)
    cur_x_vec = zeros(1,2*N);
    cur_x_vec(1:2:end) = min(g_vec(i,:),1);
    cur_x_vec(2:2:end) = cur_x_vec(1:2:end);
    cur_x_vec = repmat(cur_x_vec, num_allele_vecs(i), 1);
    hetero_inds = find(g_vec(i,:) == 1);
    hetero_inds = mat_into_vec([hetero_inds'*2-1 hetero_inds'*2]')';
    cur_x_vec(:,hetero_inds) = tmp_all_x_vec(1:num_allele_vecs(i),(end-2*log2(num_allele_vecs(i))+1):end);
    %    cur_x_vec(:,hetero_inds*2) = 1-cur_x_vec(:,hetero_inds*2-1);
    x_vec(ctr:ctr+num_allele_vecs(i)-1,:) = cur_x_vec;
    x_probs(ctr:ctr+num_allele_vecs(i)-1) = g_probs(i) / num_allele_vecs(i);
    ctr=ctr+num_allele_vecs(i);
end



