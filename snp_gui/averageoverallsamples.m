% Written by Or Zuk 5/2007
%
% This function computes the average copy number
% and also builds a matrix containing all genotypes
% and all (total) copy numbers. It does not output anything
% but saves the results into a file
%
% Input:
% user_dir - working directory with data files
% SampleNames - nameas of samples
% chip_type - type of chip
%
function AverageOverAllSamples(user_dir, SampleNames, chip_type)
num_samples = length(SampleNames);
ALL_MAT = {};

AllDispStruct = {};  AllDispStruct.SampleName = 'All Samples';
AllDispStruct.ChipName = chip_type;

for i=1:24
    AllDispStruct.Chrom{i} = []; % Start with an empty field
end
% Loop over sampels
for cur_sample=1:num_samples
    sample_name = SampleNames{cur_sample}
    display_file_name = fullfile(user_dir, 'display',[sample_name '_' chip_type '_disp']);
    load(display_file_name);

    if(cur_sample==1)
        AllDispStruct.GenomeBuild = DispStruct.GenomeBuild;
        num_snps = 0;
        for i=1:24
            if(~isempty(DispStruct.Chrom{i}))
                num_snps = num_snps + length(DispStruct.Chrom{i}.Locs);
            end
        end
        ALL_MAT.data_genotype_AB = zeros(num_samples, num_snps, 'single');
        ALL_MAT.chr_vec = zeros(1, num_snps);
        ALL_MAT.loc = zeros(1, num_snps);
    end
    
    chr_start_ind = 1; chr_end_ind=0;
    for i=1:24 % loop over chromosomes
        if(~isempty(DispStruct.Chrom{i}))
            chr_end_ind = chr_end_ind + length(DispStruct.Chrom{i}.Locs);
            if(cur_sample==1)
                AllDispStruct.Chrom{i}.Locs = DispStruct.Chrom{i}.Locs;
                AllDispStruct.Chrom{i}.Data = zeros(size(DispStruct.Chrom{i}.Data,1), size(DispStruct.Chrom{i}.Data,2));
                AllDispStruct.Chrom{i}.Genotypes = zeros(length(DispStruct.Chrom{i}.Genotypes), 1);
                ALL_MAT.chr_vec(chr_start_ind:chr_end_ind) = i;
                ALL_MAT.loc(chr_start_ind:chr_end_ind) = DispStruct.Chrom{i}.Locs';
            end
            % In genotypes we only distinguish between hetro and homo zygoce
            AllDispStruct.Chrom{i}.Genotypes = AllDispStruct.Chrom{i}.Genotypes + DispStruct.Chrom{i}.Genotypes;
            AllDispStruct.Chrom{i}.Data = AllDispStruct.Chrom{i}.Data + DispStruct.Chrom{i}.Data;

            ALL_MAT.data_genotype_AB(cur_sample,chr_start_ind:chr_end_ind) = DispStruct.Chrom{i}.Genotypes;
            ALL_MAT.data_A(cur_sample,chr_start_ind:chr_end_ind) = single(DispStruct.Chrom{i}.Data(:,1));
            ALL_MAT.data_B(cur_sample,chr_start_ind:chr_end_ind) = single(DispStruct.Chrom{i}.Data(:,2));
            chr_start_ind = chr_end_ind+1;
        end
    end
end

for i=1:24
    if(~isempty(AllDispStruct.Chrom{i}))
        AllDispStruct.Chrom{i}.Data = AllDispStruct.Chrom{i}.Data ./ num_samples; % Normalize by number of samples
        AllDispStruct.Chrom{i}.Segments = load_chr_bands(i); % Band's start&end

        % Get band SNPs average copy number
        for band=1:size(AllDispStruct.Chrom{i}.Segments, 1)

            BandInds = find( (AllDispStruct.Chrom{i}.Locs > AllDispStruct.Chrom{i}.Segments(band,1)) & ...
                (AllDispStruct.Chrom{i}.Locs <= AllDispStruct.Chrom{i}.Segments(band,2)) );
            if(~isempty(BandInds))
                AllDispStruct.Chrom{i}.Segments(band,3) = 2*mean(mean(AllDispStruct.Chrom{i}.Data(BandInds,:)));
            else
                AllDispStruct.Chrom{i}.Segments(band,3) = -1;
            end
        end
        AllDispStruct.Chrom{i}.Segments = AllDispStruct.Chrom{i}.Segments(AllDispStruct.Chrom{i}.Segments(:,3)>=0, :);
    end
end


% Create also a display file representing the average over all patients ..
ALL_MAT.SampleNames = SampleNames;
save(fullfile(user_dir, 'hmm_out', 'AllSamplesAverage.mat'), 'ALL_MAT');

DispStruct = AllDispStruct;
save(fullfile(user_dir, 'display', 'AllSamplesAverage_disp.mat'), 'DispStruct');


% The HAPMAP samples:
%% SampleNames = {'NA06985_Hind_B5_3005533affybench', 'NA07029_Hind_A9_4000092affybench',	'NA07357_Hind_B10_3005533affybench', ...
%% 'NA06991_Hind_B6_3005533affybench', 'NA07034_Hind_B1_4000092affybench', 'NA10830_Hind_E6_4000092affybench',	...
%% 'NA06993_Hind_B4_4000092affybench',	'NA07048_Hind_B3_4000092affybench',	'NA10831_Hind_E9_4000092affybench', ...
%% 'NA06994_Hind_A7_3005533affybench',	'NA07055_Hind_B2_4000092affybench',	'NA10835_Hind_E12_4000092affybench', ...
%% 'NA07000_Hind_A8_3005533affybench',	'NA07056_Hind_A11_4000092affybench',	'NA07019_Hind_A12_4000092affybench', ...
%% 'NA07345_Hind_B11_3005533affybench',	'NA07022_Hind_A10_4000092affybench',	'NA07348_Hind_B12_3005533affybench'};
%
