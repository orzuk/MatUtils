function RunManySampleFilesToOneAllSampleFile()

%% user_dir = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\CEU';
user_dir = 'F:\Or\hapmap\CEU_test'; % eytan14-pc
chip = 'xba';

AssignAllGlobalConstants;
SampleNames = HapmapSampleNames{CEU};
NumSamples = length(SampleNames);
for i=1:NumSamples
    SampleNames{i} = ['CEU_' SampleNames{i} '_n'];
end

ManySampleFilesToOneAllSampleFile(user_dir, chip, SampleNames, NumSamples);
