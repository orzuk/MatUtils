function  TwoChipsCopyChangesScore = CompareTwoChipsAbberations_script()

%user_dir = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\Leukemia';
user_dir = 'E:\public\SNP_GUI\SNP_GUI\data\Leukemia';

chip_types = {'hind', 'xba'};
chroms = [1:23];

beta = 100;
f_name_str = ['beta_' num2str(beta)];
display_dir = ['display_beta_' num2str(beta)];
chip1 = load(fullfile(user_dir, display_dir, ['all_samples_hmm_out_' chip_types{1} '.mat']), 'sample_names');
chip2 = load(fullfile(user_dir, display_dir, ['all_samples_hmm_out_' chip_types{2} '.mat']), 'sample_names');
sample_names{1} = intersect(chip1.sample_names, chip2.sample_names); 

TwoChipsCopyChangesScore{1} = CompareTwoChipsAbberationsSamples(user_dir, chip_types, display_dir, sample_names{1}, chroms, f_name_str);

beta = 10000;
f_name_str = ['beta_' num2str(beta)];
display_dir = ['display_beta_' num2str(beta)];
chip1 = load(fullfile(user_dir, display_dir, ['all_samples_hmm_out_' chip_types{1} '.mat']), 'sample_names');
chip2 = load(fullfile(user_dir, display_dir, ['all_samples_hmm_out_' chip_types{2} '.mat']), 'sample_names');
sample_names{2} = intersect(chip1.sample_names, chip2.sample_names); 
TwoChipsCopyChangesScore{2} = CompareTwoChipsAbberationsSamples(user_dir, chip_types, display_dir, sample_names{2}, chroms, f_name_str);

[sample_names_both I J] = intersect(sample_names{1}, sample_names{2});

figure; plot(TwoChipsCopyChangesScore{1}(I), TwoChipsCopyChangesScore{2}(J), '.');




