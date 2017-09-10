function ManySampleFilesToOneAllSampleFile(user_dir, chip, SampleNames, NumSamples)

AssignAllGlobalConstants;


for i=1:NumSamples
    load_file = i
    sample_str = [SampleNames{i} '_' lower(chip) '.mat'];
    load(fullfile(user_dir, sample_str), 'snp_ids', 'copy_num_vec', 'allele_ratio_vec');
    
    [NormalizedSNPsCopyMatA(:,i) NormalizedSNPsCopyMatB(:,i)] = RatioToCopyMats(copy_num_vec, allele_ratio_vec);
end

% Addition: Save also one big matrix with everything. Save A and B intensities
save(fullfile(user_dir, ['TempAllSamplesMat_' lower(chip) '.mat']), 'snp_ids', 'SampleNames', 'chip', ...
    'NormalizedSNPsCopyMatA', 'NormalizedSNPsCopyMatB');



