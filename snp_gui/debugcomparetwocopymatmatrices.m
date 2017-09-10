function DebugCompareTwoCopymatMatrices(ALL_MAT)

AssignAllGlobalConstants;


for hapmap_population = 1:3
    ALL=load(fullfile('F:\Or\hapmap\', [pop_str_vec{hapmap_population} '_test'], 'OldAllSamplesMat_xba.mat'));
    DOUBLE_ALL=load(fullfile('F:\Or\hapmap\', [pop_str_vec{hapmap_population} '_test'], 'AllSamplesMat_xba.mat'));

    ONE_SAMPLE = load(fullfile('F:\Or\hapmap\', [pop_str_vec{hapmap_population} '_test'], 'CEU_NA07029_n_xba.mat'));
    [S I J] = intersect(ALL.snp_id_xba, ALL_MAT.snp_ids);
        
    min(min(ALL.NormalizedSNPsCopyMatA(I,:) - ALL_MAT.data_A((hapmap_population-1)*90+1:hapmap_population*90,J)')).^2 
end

'F:\Or\hapmap\CEU_test\SingleAllSamplesMat_xba.mat', 