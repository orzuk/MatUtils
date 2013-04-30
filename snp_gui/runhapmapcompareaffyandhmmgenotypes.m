%% comp_data = 'AffyBenchmark_data';
%comp_data = 'Luekemia';

AssignAllGlobalConstants();

comp_data = 'AffyBenchmark_data';  % This is the hapmap data
chip_type = 'hind';
chroms = [1:23];
num_samples = 2;

hapmap_dir = '\\eytan14-pc\Or\hapmap\Genotype_data';
hapmap_population = CEU; hapmap_version = 'r22';


if(strcmp(comp_data,'AffyBenchmark_data'))
    user_dir = fullfile('\\eytan14-pc\Or\hapmap', [pop_str_vec{hapmap_population} '_test']);   % E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\';

    comp_hap_affy = 0; comp_hap_hmm = 1; comp_affy_hmm = 0;% Which ones to compare ?
    use_strand_vec = [0 1 0]; % affy and hmm are different strands
else
    user_dir = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\Leukemia\';
    % Try two sample together
    comp_hap_affy = 0; comp_hap_hmm = 0; comp_affy_hmm = 1;% Which ones to compare ?
    use_strand_vec = [0 1 1]; % affy and hmm are different strands
end

% Old version:
% ErrorStruct = HapMapCompareAffyAndHMMGenotypes(FirstSample, chip_type, user_dir, ...
%     comp_hap_affy, comp_hap_hmm, comp_affy_hmm, use_strand_vec);

ErrorStruct = HapMapCompareAffyAndHMMGenotypes(num_samples, chip_type, user_dir, chroms, ...
    comp_hap_affy, comp_hap_hmm, comp_affy_hmm, use_strand_vec, ...
    hapmap_dir, hapmap_population, hapmap_version, RLMM_population)


%% save E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\GenotypeErrorStruct.mat ErrorStruct;

