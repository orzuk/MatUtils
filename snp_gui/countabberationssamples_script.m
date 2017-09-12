function  CountAbberationsSamples_script()

user_dir = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\Leukemia';

chip_type = 'xba';
load(fullfile(user_dir, 'display', ['all_samples_hmm_out_' chip_type '.mat']), 'sample_names');
SampleNames = sample_names;
f_name_str = 'beta_100';
chroms = [1:23];
[NumCopyChanges CopyNumCountsMat] = CountAbberationsSamples(user_dir, SampleNames, chip_type, chroms);

save(fullfile(user_dir, 'display', ['count_aberrations_' chip_type f_name_str '.mat']), 'NumCopyChanges', 'CopyNumCountsMat', 'SampleNames');
