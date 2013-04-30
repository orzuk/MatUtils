% Compare paired of sampled - tumor and normal, and see how many of the
% genotypes are the same (according to Affymetrix readings)
% The inputs are:
% user_dir - directory with the file 'AllSamplesAverage.mat'
% SampleNames - The sample names 
% 
% Output: 
% ErrorStruct - contains the number of SNPs with different genotypes
% between Normal and Tumor.
function ErrorStruct = CompareTumorAndNormalGenotypes(user_dir, SampleNames, chip_type)

% First process the sample names to find the matches 
TempSampleNames = SampleNames;

for i=1:length(SampleNames)
    TempSampleNames{i} = SampleNames{i}(1:end-1);
end
[CoupleNames I J] = unique(TempSampleNames);

NumCouples = length(CoupleNames);

% loop over the couples 
ErrorStruct.AffyDiffSNPs = zeros(1,NumCouples); ErrorStruct.HMMDiffSNPs = zeros(1,NumCouples);
for i=1:NumCouples
    do_hmm_i=i
    Normal = load(fullfile(user_dir, 'display', [CoupleNames{i}, 'n_' chip_type '_disp.mat'])); % 'n' for normal
    Tumor = load(fullfile(user_dir, 'display', [CoupleNames{i}, 'd_' chip_type '_disp.mat'])); % 'd' for disease

    ErrorStruct.HMMDiffSNPs(i) = 0;
    for j=1:24
        if(~isempty(Normal.DispStruct.Chrom{j}))
            TEMP_NORMAL = bitget(Normal.DispStruct.Chrom{j}.Genotypes,1) + bitget(Normal.DispStruct.Chrom{j}.Genotypes,2);
            TEMP_TUMOR = bitget(Tumor.DispStruct.Chrom{j}.Genotypes,1) + bitget(Tumor.DispStruct.Chrom{j}.Genotypes,2);
           
            ErrorStruct.HMMDiffSNPs(i) = ErrorStruct.HMMDiffSNPs(i) + sum(TEMP_NORMAL ~= TEMP_TUMOR); % here AB and BA are the same ...
%%                sum(Normal.DispStruct.Chrom{j}.Genotypes ~= Tumor.DispStruct.Chrom{j}.Genotypes);
        end
    end
end

get_affy_genotypes = 0;  %take the genotypes from the hmm and not from affymetrix.
for i=1:NumCouples
    do_i = i
    ErrorStruct.CoupleNames{i} = CoupleNames{i};
    if(get_affy_genotypes == 1) % here takes the original genotypes from affymetrix
        Normal = load(fullfile(user_dir, [CoupleNames{i}, 'n_' chip_type '.mat'])); % 'n' for normal
        Tumor = load(fullfile(user_dir, [CoupleNames{i}, 'd_' chip_type '.mat'])); % 'd' for disease

        ErrorStruct.SampleNames{i} = CoupleNames{i};
        ErrorStruct.AffyDiffSNPs(i) = sum(Normal.genotype_vec ~= Tumor.genotype_vec);    

    else % get genotypes from hmm output - this can be actually ignored if we don't have affymetrix genotypes ... 
        Normal = load(fullfile(user_dir, 'display', [CoupleNames{i}, 'n_' chip_type '_disp.mat'])); % 'n' for normal
        Tumor = load(fullfile(user_dir, 'display', [CoupleNames{i}, 'd_' chip_type '_disp.mat'])); % 'd' for disease
        
        Normal.genotype_vec = []; Tumor.genotype_vec = [];
        
        for j=1:length(Normal.DispStruct.Chrom)
            if(~isempty(Normal.DispStruct.Chrom{j}))
                Normal.genotype_vec = [Normal.genotype_vec Normal.DispStruct.Chrom{j}.Genotypes'];
                Tumor.genotype_vec = [Tumor.genotype_vec Tumor.DispStruct.Chrom{j}.Genotypes'];
            end
        end
    end
    if(i == 1)
        ErrorStruct.NumSNPs = length(Normal.genotype_vec);
    end
end





