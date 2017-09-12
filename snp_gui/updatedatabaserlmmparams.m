% Written by Or Zuk 7/2007
%
% A function for updating the Gaussian moments parameters to be saved in the
% database for the genotyping algorithm. The method is simply using a training
% set and computing the moments based on it as described in the RLMM paper.
%
% The inputs:
% user_dir - working directory with input files
% SampleNames - names of samples in the training set
% LDStruct - structure containing Linkeage-Disequilibrium statistics from Hapmap (not needed ?)
% SNPChipAnnotStruct - Annotations for the specific SNP chip
% HMMParamsStruct - parameters for the HMM running (Needed ?)
% hapmap_population - which hapmap population do we want to use in order to
% determine the parameters (CEU, JPT_CHB, YRI or all of them together ?)
% hapmap_version - which hapmap_version do we use for reading genotypes?
% MinSparse  - Mimimal number of SNPs with a given genotype to use for moment estimation
%
% The outputs:
% OutputFilesNames - the name of the files in which the function saves the Gaussian Parameters.
function OutputFilesName = ...
    UpdateDatabaseRLMMParams(user_dir, SampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct, ...
    hapmap_dir, hapmap_population, hapmap_version, MinSparse, PsuedoCount)

% Various data types
AssignAllGlobalConstants();
HMMParamsStruct.use_affy_genotypes = 0; % Should be eliminated - input as a button ...
chip_type = lower(SNPChipAnnotStruct.chip);
num_samples = length(SampleNames);
chroms = 1:24;

ALL_MAT = {}; TEMP_DEBUG_FLAG=0;  % Try not to use this anymore
%%if(TEMP_DEBUG_FLAG)
%%    load('TEMP_ALL_CEU_HIND.mat'); % Very specific - remove in the end


% get the samples genders
load(fullfile('..', 'database', 'AllHapmapGenders.mat'));

for i=1:length(AllHapMapSamples.Name)
     AllHapMapSamples.Name{i} = [ AllHapMapSamples.Name{i} '_n'];
end

ctr=0;
if(hapmap_population == ALL_POPS) % Take everything together - Currently NOT WORKING !!!
    for hapmap_populations2 = [CEU JPT_CHB YRI]
        doing_population = pop_str_vec{hapmap_populations2}
        ttt = cputime;

        CurSampleNames = SampleNames(strmatch(pop_str_vec{hapmap_populations2}, SampleNames));
        cur_num_samples = length(CurSampleNames);
        [cur_samples_inter I J] = intersect(CurSampleNames, AllHapMapSamples.Name);
        CurSamplesGenders(I) = AllHapMapSamples.Gender(J);
        SamplesGenders(ctr+[1:cur_num_samples]) = CurSamplesGenders;
        %        [A_COPY B_COPY AB_GENOTYPES] =UpdateOnePopulation(hapmap_populations2, fullfile(user_dir, [pop_str_vec{hapmap_populations2} '_test']), ...
        %            CurSampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct, hapmap_dir, hapmap_version);

        [ALL_MAT.data_A(ctr+[1:cur_num_samples],:) ALL_MAT.data_B(ctr+[1:cur_num_samples],:) ...
            ALL_MAT.genotypes_AB(ctr+[1:cur_num_samples],:) ALL_MAT.snp_ids ALL_MAT.strand ALL_MAT.chr]  = ...
            UpdateOnePopulation(hapmap_populations2, fullfile(user_dir, [pop_str_vec{hapmap_populations2} '_test']), chip_type, ...
            CurSampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct, hapmap_dir, hapmap_version);
        ctr=ctr+cur_num_samples;
        time_passed = cputime-ttt

    end
else
    [cur_samples_inter I J] = intersect(SampleNames, AllHapMapSamples.Name);
    % Deal with special sample with a bad name for the JPT_CHB population
    BAD_SPECIAL_SAMPLE = strmatch('JPT_CHB_NA18996_n', SampleNames);
    if(~isempty(BAD_SPECIAL_SAMPLE))
        BAD_SPECIAL_SAMPLE_HAPMAP = strmatch('JPT_CHB_NA19012_n', AllHapMapSamples.Name);
        if(~isempty(BAD_SPECIAL_SAMPLE_HAPMAP))
            I = [I BAD_SPECIAL_SAMPLE]; J = [J BAD_SPECIAL_SAMPLE_HAPMAP];
        end
    end
    SamplesGenders(I) = AllHapMapSamples.Gender(J);
    [ALL_MAT.data_A ALL_MAT.data_B ALL_MAT.genotypes_AB ALL_MAT.snp_ids ALL_MAT.strand ALL_MAT.chr] = ...
        UpdateOnePopulation(hapmap_population, user_dir, chip_type, SampleNames, ...
        LDStruct, SNPChipAnnotStruct, HMMParamsStruct, hapmap_dir, hapmap_version);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Gaussian Parameters for  RLMM:                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RLMM = RLMM_LearnGaussianParams(ALL_MAT.data_A', ALL_MAT.data_B', ALL_MAT.genotypes_AB', ALL_MAT.strand, ALL_MAT.chr, SamplesGenders, MinSparse, PsuedoCount);
RLMM.snp_ids = ALL_MAT.snp_ids;

OutputFilesName = fullfile('../database', ['RLMM_' pop_str_vec{hapmap_population} '_' chip_type '_' genome_assembly '.mat'])
save(OutputFilesName, 'RLMM');

%% TestRLMMAdjust;













