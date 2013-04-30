%function  TwoChipsCopyChangesScore = CompareTwoChipsAbberationsSamples(user_dir, chip_types, display_dir, sample_names, chroms, f_name_str)
function  TwoChipsCopyChangesScore = CompareTwoChipsAbberationsSamples(user_dir, chip_types, display_dir, sample_names, chroms, f_name_str)


num_samples = length(sample_names)
TwoChipsCopyChangesScore = zeros(1, num_samples);
for s = 1:num_samples
    sample_name = sample_names{s}
    for c = 1:length(chip_types)
        chip_type = chip_types{c};
        display_file_name = fullfile(display_dir, [sample_name '_' chip_type '_disp']);
        D{c} = load(fullfile(user_dir, display_file_name),'DispStruct'); % load the display structure
    end
    TwoChipsCopyChangesScore(s) = CompareTwoDispStructsAbberations(D{1}.DispStruct, D{2}.DispStruct, chroms);
end


save(fullfile(user_dir, display_dir, ['comp_copy_jump_aberr_' chip_types{1} '_' chip_types{2} '_' f_name_str '.mat']), ...
    'TwoChipsCopyChangesScore', 'sample_names');
