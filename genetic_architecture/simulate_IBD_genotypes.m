% Simulate genotypes for a set of individuals with a given IBD-sharing matrix
%
% Input:
% IBD_mat - matrix of IBD sharing values (between zero and one) 
% f_vec - vector of allele frequencies
% iters - how many instances to simulate
%
% Output:
% SNP_mat - matrix of SNPs-by-individuals
% corr_mat - matrix of pairwise correlations between individuals based on SNPs
%
function [SNP_mat corr_mat] = ...
    simulate_IBD_genotypes(IBD_mat, f_vec, iters) % simulate genotype vectors for entire family

use_matlab = 0; % 0 - compute conversion to Gaussian in matlab, 1 - use sampleDichGauss01 (faster) 
n = length(IBD_mat); % number of individuals
num_snps = length(f_vec); % number of SNPs 
SNP_mat = zeros(num_snps,n);

[sorted_f_vec sort_perm] = sort(f_vec);
new_f_inds = [1 1+find(sorted_f_vec(1:end-1) ~= sorted_f_vec(2:end))' num_snps+1]; % find change in f AFTER sorting
num_to_simulate_vec = diff(new_f_inds);
num_unique_snps = length(new_f_inds)-1;
new_ctr=0;

for k=1:num_unique_snps % loop on all uniques SNPs
    if(mod(k,50)==0)
        simulate_snp = k
        total_snps = num_snps
    end
    if(use_matlab)
        [corr_mat x_f_vec] = IBD_mat_to_corr_mat(IBD_mat, f_vec(k)); % take the first allele frequency
        z = simulate_correlated_gaussians(iters, length(corr_mat), corr_mat);
        SNP_mat(k,:) = z<x_f_vec;
    else
        
        % % %         if((k>1) && (sorted_f_vec(k) == sorted_f_vec(k-1))) % already have the moments
        % % %             if(mod(k,50)==0)
        % % %                 simulate_same_f = k
        % % %             end
        % % %                     SNP_mat(k,:) = sampleDichGauss01(mu, sigma,iters, 1); % ,already_computed,acc)
        % % %            % z = simulate_correlated_gaussians(iters, n, sigma);
        % % %            % SNP_mat(k,:) = z<mu(1);
        % % %
        % % %         else  % compute Gaussian moments (new inds)
        simulate_new_f = k
        binary_corr_mat = sorted_f_vec(new_f_inds(k)).*(1-sorted_f_vec(new_f_inds(k))).*IBD_mat;
        [SNP_mat(new_f_inds(k),:), mu, sigma] = ...
            sampleDichGauss01(repmat(sorted_f_vec(new_f_inds(k)), 1, n), binary_corr_mat,iters); % ,already_computed,acc)
        if(~isposdef(sigma))
            error('Error! Matrix not pos-def !!!\n');
        end
        new_ctr=new_ctr+1; % update new counter

% %         for jj=new_f_inds(k)+1:new_f_inds(k+1)-1
% %             SNP_mat(jj,:) = ...
% %                 sampleDichGauss01(mu, sigma, 1, 1); % ,already_computed,acc)
% %         end
        SNP_mat(new_f_inds(k)+1:new_f_inds(k+1)-1,:) = ...
            sampleDichGauss01(mu, sigma,num_to_simulate_vec(k)-1, 1)'; % ,already_computed,acc)
        
        % % %         end % if new f_index
    end
end % loop on unique SNPs 

if(~use_matlab)
    inv_sort_perm = inv_perm(sort_perm);
    SNP_mat = SNP_mat(inv_sort_perm,:); % Sort back to original order
    corr_mat = sigma; % take
end



