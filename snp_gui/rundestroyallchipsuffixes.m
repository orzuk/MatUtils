% Run the short utility that tries to eliminate all the '_hind' created by older
% version of the program WITHOUT running everything again


destroy_dir{1} = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\CEU';
destroy_dir{2} = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\YRI';
destroy_dir{3} = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\JPT_CHB';
destroy_dir{4} = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\Leukemia'

for j=1:length(destroy_dir) 
    j_is = j
    j_is = j
    j_is = j
    j_is = j
    j_is = j
    j_is = j
    j_is = j
    j_is = j
    j_is = j
    SpecialInds = DestroyAllChipSuffixes(destroy_dir{j})
end
