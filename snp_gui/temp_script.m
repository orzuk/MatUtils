AssignAllGlobalConstants;
load('F:\Or\LeukemiaDchip\HD78_9_n_xba.mat')

[CopyMatA CopyMatB] = RatioToCopyMats(copy_num_vec_xba, allele_ratio_vec_xba); 

hapmap_population = CEU; chip_type = 'xba'; PsuedoCount=40;
RLMM_Correction = RLMM_AdjustParamsToCurrentData(CopyMatA, CopyMatB, snp_id_xba, hapmap_population, chip_type, PsuedoCount)
