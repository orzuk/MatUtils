% Calculate the frequency of every SNP in a given population. We also
% calculate the pairwise statistics (Linkage Disequilibrium), and also
% return the snp ids vector
function [A_freq_vec, pairwise_freq_mat, data_snp_ids] = calc_allele_freq(user_dir, chip_type)

load(fullfile(user_dir, 'display', ['all_samples_hmm_out_' chip_type '.mat']), 'genotype_mat', 'data_snp_ids');

% Change AB/BA to 1 and BB to 2 
genotype_mat(genotype_mat == 2) = 1;
genotype_mat(genotype_mat == 3) = 2;

A_freq_vec = single(1 - mean(genotype_mat, 2) ./ 2); % get the singleton frequencies

% Now get the pairwise frequencies
pairwise_freq_mat = zeros(4,length(A_freq_vec));

%[chr, chr_locs, pres_ind] = get_chrom_locs(data_snp_ids, chip_type); % get the locations

% sort the SNPs according to the chromosomes, and then the locations
%[sorted_locs I J] = sort(chr_locs);


