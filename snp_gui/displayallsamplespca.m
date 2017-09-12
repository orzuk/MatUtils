% Read the files containing all the stuff on samples and display it.
% We display either accoding to genotypes, OR according to total copy
% number (A+B). For the labels we take the last letter, 
% so we can distinguish for example between disease and healthy,
% or between different populations
% The inputs are:
% user_dir - directory with the file 'AllSamplesAverage.mat'
% data_flag - 0: copy number, 1 - genotypes
% data_labels - some labels (e.g.  samples names)
% chrom - which chromosome to show (-1 is all chromsomes)
function DisplayAllSamplesPCA(user_dir, chip_type, data_flag, data_labels, chrom, NumRandSNPs)

% load the data
load( fullfile(user_dir, 'hmm_out', 'AllSamplesAverage.mat') );

% extract the last letter - this is the (old !!) label setting the color in the PCA plot !! 
data_labels_last = data_labels;
for i=1:length(data_labels)
    data_labels_last{i} = data_labels{i}(end);
end

if(data_flag == 0) % Show copy number (we take raw data - why not HMM output?)
    ALL_MAT = load(fullfile(user_dir, ['AllSamplesMat_' chip_type '.mat']), 'snp_ids', 'NormalizedSNPsCopyMatA', 'NormalizedSNPsCopyMatB');
    ALL_MAT.data_A = ALL_MAT.NormalizedSNPsCopyMatA'; ALL_MAT = rmfield(ALL_MAT, 'NormalizedSNPsCopyMatA');
    ALL_MAT.data_B = ALL_MAT.NormalizedSNPsCopyMatB'; ALL_MAT = rmfield(ALL_MAT, 'NormalizedSNPsCopyMatB'); 
else % show genotypes
    ALL_MAT = load(fullfile(user_dir, 'display', ['all_samples_hmm_out_' chip_type '.mat']), ...
        'data_snp_ids', 'genotype_mat', 'sample_names');
    ALL_MAT.snp_ids = ALL_MAT.data_snp_ids; ALL_MAT = rmfield(ALL_MAT, 'data_snp_ids');
    ALL_MAT.data_genotype_AB = ALL_MAT.genotype_mat'; ALL_MAT = rmfield(ALL_MAT, 'genotype_mat');
end
ALL_MAT.chip_type = chip_type;

% call PCA display function
% We assume here that the ALL_MAT structure contains the following fields:
% snp_ids (obvious)
% data_a, data_b - for the case that data_flag is copy #
% data_genotype_AB - for the case where data_flag is genotypes 
DisplaySamplesPCA(ALL_MAT, data_flag, data_labels_last, chrom, NumRandSNPs); 




