AssignAllGlobalConstants;
chip_type = 'xba';

CurSample = load('F:\Or\Leukemia\HD82_3_d_xba.mat'); % load new - our normalization % Eytan14-pc


load('F:\Or\Leukemia\HD78_9_n_xba.mat'); % load new - our normalization % Eytan14-pc
%% load('E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\Leukemia\HD78_9_n_xba.mat'); % load new - our normalization

allele_ratio_vec = max(EPSILON, min(allele_ratio_vec, 99999999999));
ret_strand = load_snp_ids_strand(snp_ids, chip_type);

BadInds = find(isnan(allele_ratio_vec));
GoodInds = setdiff(1:length(copy_num_vec),BadInds);
[CopyMatA CopyMatB] = RatioToCopyMats(copy_num_vec(GoodInds), allele_ratio_vec(GoodInds)); 

hapmap_population = ALL_POPS;  PsuedoCount=40;
ret_strand = ret_strand(GoodInds);
minus_strand = strmatch('-', ret_strand);
[CopyMatA(minus_strand, :) CopyMatB(minus_strand, :)] = swap(CopyMatA(minus_strand, :), CopyMatB(minus_strand, :));

RLMM_Correction = RLMM_AdjustParamsToCurrentData(CopyMatA, CopyMatB, snp_ids(GoodInds), ...
    hapmap_population, chip_type, PsuedoCount)
