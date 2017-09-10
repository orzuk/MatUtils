% Here we implement the RLMM package of Rabbee&Speed (Bioinformatics 2006),
% which was originally written in R, in Matlab.
% This file does the raw normalization. There are two steps:
% 1. Quantile normalization.
% 2. Robust Linear fit modelling.
%
% The input to the function is a bunch of .cell files lying in the same
% directory. (Note that in the original RLMM they require .raw files which
% are .cell files after some transformation.). There is no output variable
% to the function - it saves the normalized samples into output files
% The input:
%
% SampleNames - Name of samples.
% NormalsInds - The indexes of normal samples (used for normalization)
% FemaleInds - The indexes of females (for chrom. X normalization)
% arrays - what the hell is this ?
% chip - chip name
% load_cel_files_flag - 1 if we load the .cel files (this must be done in the first time,
% but it takes lots of time so its better not to do it afterwards).
% prep_block_flag - 1 if we want to prepare the blocks (also takes time),
% if it is 0 we use the prepared block files
%
function [FemaleInds] =RLMM_Normalize(user_dir, SampleNames, NormalsInds, FemaleInds, arrays, chip, ...
    load_cel_files_flag, prep_block_flag, gender_not_given_vec)

AssignAllGlobalConstants();

std_outlier = 3; % for calculating the % array outlier
outlier_thresh = 0.1;
remove_outlier_from_normals = 1;
%remove_outlier_from_normals = 0;
%std_iterate_flag = 1;
std_iterate_flag = 0;
quantile_norm_flag = 0;
%quantile_norm_flag = 1;
median_std_flag = 1;
%median_std_flag = 0;

% % Load Chip Annotations - We need the X,Y chromosmes from here.
SNPChipAnnotStruct = load(fullfile('..\database', [chip '_annot_data_' genome_assembly '.mat']), 'snp_ids', 'chr_vec', 'chr_loc_vec');


%% load_cel_files_flag = 1; load_cel_files_flag = 0; prep_block_flag = 1;
if(load_cel_files_flag==1)
    prep_block_flag = 1;
end

NumSamples = length(SampleNames);
cdf_f_name = get_cdf_f_name(chip);
CDFStruct=affyread(cdf_f_name);
probesets = get_probesets(chip);

median_intensity_mat = zeros(NumSamples, 2); % keeps the median intensity before and after quantile normalization

if(~load_cel_files_flag)
    % load snp_ids in case we don't load CEL files this run
    sample_str = [SampleNames{1} '_' lower(chip) '.mat'];
    %    load(fullfile(user_dir, sample_str), snp_id_str);
    load(fullfile(user_dir, ['AllSamplesMat_' lower(chip) '.mat']), 'snp_ids');
    NumSNPs = length(snp_ids);
end

if(strcmp(lower(chip), 'hind') | strcmp(lower(chip), 'xba'))
    pa_vec = [1:10]; ma_vec = [21:30]; pb_vec = [11:20]; mb_vec = [31:40];
else
    pa_vec = [1:14]; ma_vec = [29:42]; pb_vec = [15:28]; mb_vec = [43:56];
end
num_pm = length(pa_vec);
pm_ind = [1:num_pm*2];
num_probes = 4*num_pm;

if(load_cel_files_flag)
    % 1. Quantile normalization.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop and read all files
    for i=1:NumSamples
        sample_name = SampleNames{i};
        cel_data_f = fullfile(user_dir, [sample_name, '_CEL_data_' lower(chip) '.mat']);
        % check if this cel file was read before
        if ~exist(cel_data_f,'file')

            if strcmpi(chip, 'sty')
                clear probes snp_ids data
                save tmp_RLMM_Normalize.mat
                clear
                load tmp_RLMM_Normalize.mat arrays i CDFStruct probesets chip
            end
            [probes snp_ids data]=load_cel_probes(arrays{i}, CDFStruct, probesets);
            if strcmpi(chip, 'sty')
                load tmp_RLMM_Normalize.mat
                delete tmp_RLMM_Normalize.mat
            end

            % parse the probes variable
            ProbesPerm = GetProbesPerm(probes);
            data = data(:,ProbesPerm);
            if(strcmp(lower(chip), 'hind') | strcmp(lower(chip), 'xba'))
                if(size(data, 2) ~= num_probes)
                    data_vec = mat_into_vec(data');
                    not_nan_ind = find(1-isnan(data_vec(:))); % Get 40 (num_probes) not-Nan out of 56
                    data = vec_into_mat(data_vec(not_nan_ind), num_probes)'; %
                end
            end
            save( cel_data_f, 'data', 'snp_ids');
        else
            if(i==1)
                load( cel_data_f, 'data', 'snp_ids');
            else
                load( cel_data_f, 'data');
            end
        end


        if(i == 1)
            NumSNPs = size(data,1);
            allele_ratio_mat = zeros(NumSNPs, NumSamples, 'single');
        end
        
        allele_ratio_mat(:,i) = calc_probe_allele_ratio(data, 1, pa_vec, ma_vec, pb_vec, mb_vec);
        median_intensity_mat(i,1) = median(mat_into_vec(data(:, pm_ind)));

        % Sort to get quantiles
        [R_sorted R_inds] = sort(mat_into_vec(data(:, pm_ind)));

        if(i == 1)
            Q = R_sorted;
        else
            Q = Q + R_sorted;
        end
        t=[sample_name ' was read' ]
    end
end

clear CDFStruct data_vec not_nan_ind 
Q = Q ./ NumSamples; % Take the mean

if(quantile_norm_flag)
    % loop again and read CEL data files
    for i=1:NumSamples
        sample_name = SampleNames{i};

        cel_data_f = fullfile(user_dir, [sample_name, '_CEL_data_' lower(chip) '.mat']);
        load(cel_data_f, 'data');
        data_pm = subtract_mm_from_pm(data, pa_vec, ma_vec, pb_vec, mb_vec);
        % Sort to get quantiles
        %        data_pm = data(:, pm_ind);
        [R_sorted R_inds] = sort(data_pm(:));
        R_sorted(R_inds) = Q;
        data_pm = vec_into_mat(R_sorted, NumSNPs);
        save( fullfile(user_dir,[sample_name, '_PM_norm1_' lower(chip) '.mat']), 'data_pm');
    end

elseif(median_std_flag)
    end_median = 2000;
    end_std = 2000;
    for i=1:NumSamples
        sample_name = SampleNames{i};

        cel_data_f = fullfile(user_dir, [sample_name, '_CEL_data_' lower(chip) '.mat']);
        load(cel_data_f, 'data');
        % normalize to have equal median and std
        %            data_pm = data(:, pm_ind);
        data_pm = subtract_mm_from_pm(data, pa_vec, ma_vec, pb_vec, mb_vec);
        %        data_pm = data;
        std_data = std(data_pm(:))
        mult_fact = end_std/std_data;
        data_pm = data_pm*mult_fact;
        median_data = median(data_pm(:))
        move_fact = end_median-median_data;
        data_pm = data_pm+move_fact;
        new_median_data = median(data_pm(:))
        new_std_data = std(data_pm(:))

        %        data_pm = data_pm(:, pm_ind);
        save( fullfile(user_dir,[sample_name, '_PM_norm1_' lower(chip) '.mat']), 'data_pm');

    end
else
    load(fullfile(user_dir, ['AlleleRatioMat_' lower(chip) '.mat'])); % load the already made allele_ratio_mat
end


% Save the allele ratio matrix - so that we don't want to need to load the .cell files again
% Sort SNPs by chromosome and chromosomal location
[Ids I J] = intersect(snp_ids, SNPChipAnnotStruct.snp_ids);
[J_sorted J_perm] = sort(J); I = I(J_perm);
snp_ids = snp_ids(I);

if(load_cel_files_flag)
    allele_ratio_mat = allele_ratio_mat(I,:);
    save(fullfile(user_dir, ['AlleleRatioMat_' lower(chip) '.mat']), 'allele_ratio_mat', 'SampleNames');
end



if(length(gender_not_given_vec)>0)
    % find chrx snps
    chrx_ind = SNPChipAnnotStruct.chr_vec(J_sorted) == 23;
    [gender_cell, female_inds] = get_gender_from_allele_ratio...
        (allele_ratio_mat(chrx_ind, gender_not_given_vec) , chip);
    if(size(FemaleInds,1)~= 1) FemaleInds = FemaleInds'; end
    if(size(gender_not_given_vec,1)~= 1) gender_not_given_vec = gender_not_given_vec'; end
    FemaleInds = [FemaleInds gender_not_given_vec(female_inds)];

end

% 2. Robust Linear fit modelling.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNPsBlockSize = min(10000, NumSNPs);
% SamplesBlockSize = min(5, NumSamples);
%SNPsBlockSize = min(2500, NumSNPs);
%SamplesBlockSize = min(20, NumSamples);
SamplesBlockSize = min(2, NumSamples);

NumSNPsBlocks = ceil(NumSNPs/SNPsBlockSize);
NumSamplesBlocks = ceil(NumSamples/SamplesBlockSize);

NormalizedSNPsMat = zeros(NumSNPs, NumSamples, 'single');
%temp_check_allele_ratio = 1;
temp_check_allele_ratio = 0;
if(temp_check_allele_ratio)
    RMA_A = zeros(NumSNPs, NumSamples, 'single');
    RMA_B = zeros(NumSNPs, NumSamples, 'single');
end

if(prep_block_flag)
    % Generate block tables for sets of SNPs and samples
    for b=1:NumSamplesBlocks
        sample_start_ind = (b-1)*SamplesBlockSize+1;
        sample_end_ind = min(b*SamplesBlockSize, NumSamples);
        num_samples_loop = sample_end_ind-sample_start_ind+1;
        big_snp_matPA = zeros(NumSNPs*num_pm, num_samples_loop, 'single');
        big_snp_matPB = zeros(NumSNPs*num_pm, num_samples_loop, 'single');
        b
        for i=sample_start_ind:sample_end_ind
            sample_name = SampleNames{i};
            load( fullfile(user_dir, [sample_name, '_PM_norm1_' lower(chip) '.mat']), 'data_pm');
            big_snp_matPA(:,i-sample_start_ind+1) = mat_into_vec(data_pm(:,1:num_pm)'); % PA
            big_snp_matPB(:,i-sample_start_ind+1) = mat_into_vec(data_pm(:,num_pm+1:num_pm*2)'); % PB
        end

        for s=1:NumSNPsBlocks
            snp_start_ind = (s-1)*SNPsBlockSize+1;
            snp_end_ind = min(s*SNPsBlockSize, NumSNPs);
            SnpsInds = [(snp_start_ind-1)*num_pm+1:snp_end_ind*num_pm];
            snp_mat_pa = big_snp_matPA(SnpsInds, :);
            snp_mat_pb = big_snp_matPB(SnpsInds, :);

            file_name = fullfile(user_dir, ['Blocks_Samples_' num2str(sample_start_ind) '_' num2str(sample_end_ind) '_SNPs_' ...
                num2str(snp_start_ind) '_' num2str(snp_end_ind) '_' lower(chip) '.mat']);
            save(file_name, 'snp_mat_pa', 'snp_mat_pb', 'sample_start_ind', 'sample_end_ind', 'snp_start_ind', 'snp_end_ind');
        end
    end
end
t=1; %#ok<NASGU>
%
%
% loop over the SNPs. For each set of SNPs, need to load all samples.
for s=1:NumSNPsBlocks
    s
    snp_start_ind = (s-1)*SNPsBlockSize+1;
    snp_end_ind = min(s*SNPsBlockSize, NumSNPs);
    num_snps_loop = snp_end_ind-snp_start_ind+1;
    big_snp_matPA = zeros(num_snps_loop*num_pm, NumSamples, 'single');
    big_snp_matPB = zeros(num_snps_loop*num_pm, NumSamples, 'single');

    for b=1:NumSamplesBlocks
        sample_start_ind = (b-1)*SamplesBlockSize+1;
        sample_end_ind = min(b*SamplesBlockSize, NumSamples);
        num_samples_loop = sample_end_ind-sample_start_ind+1;

        file_name = fullfile(user_dir, ['Blocks_Samples_' num2str(sample_start_ind) '_' num2str(sample_end_ind) '_SNPs_' ...
            num2str(snp_start_ind) '_' num2str(snp_end_ind) '_' lower(chip) '.mat']);
        load(file_name, 'snp_mat_pa', 'snp_mat_pb', 'sample_start_ind', 'sample_end_ind', 'snp_start_ind', 'snp_end_ind');

        big_snp_matPA(:, sample_start_ind:sample_end_ind) = snp_mat_pa;
        big_snp_matPB(:, sample_start_ind:sample_end_ind) = snp_mat_pb;

    end

    IndsVec = repmat([0:num_pm-1]', num_snps_loop,1);
    if(temp_check_allele_ratio)
        RMA_A(snp_start_ind:snp_end_ind,:) = single(2.^rmasummary(IndsVec ,big_snp_matPA));
        RMA_B(snp_start_ind:snp_end_ind,:) = single(2.^rmasummary(IndsVec ,big_snp_matPB));
        NormalizedSNPsMat(snp_start_ind:snp_end_ind,:) = RMA_A(snp_start_ind:snp_end_ind,:)+RMA_B(snp_start_ind:snp_end_ind,:);
    else
        NormalizedSNPsMat(snp_start_ind:snp_end_ind,:) = ...
            single(2.^rmasummary(IndsVec ,big_snp_matPA) + 2.^rmasummary(IndsVec ,big_snp_matPB));
    end

end
if(temp_check_allele_ratio)
    file_name = fullfile(user_dir, ['RMA_A_B_' lower(chip) '.mat']);
    save(file_name, 'RMA_A', 'RMA_B');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete files that are not needed anymore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:NumSamples
    sample_name = SampleNames{i};
    norm1_f_name = fullfile(user_dir,[sample_name, '_PM_norm1_' lower(chip) '.mat']);
    if(exist(norm1_f_name,'file'))
        delete(norm1_f_name);
    end
end
for s=1:NumSNPsBlocks
    snp_start_ind = (s-1)*SNPsBlockSize+1;
    snp_end_ind = min(s*SNPsBlockSize, NumSNPs);
    for b=1:NumSamplesBlocks
        sample_start_ind = (b-1)*SamplesBlockSize+1;
        sample_end_ind = min(b*SamplesBlockSize, NumSamples);
        file_name = fullfile(user_dir, ['Blocks_Samples_' num2str(sample_start_ind) '_' num2str(sample_end_ind) '_SNPs_' ...
            num2str(snp_start_ind) '_' num2str(snp_end_ind) '_' lower(chip) '.mat']);
        if(exist(file_name,'file'))
            delete(file_name);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete files that are not needed anymore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort SNPs by chromosome and chromosomal location (check that I wasn't changed).
NormalizedSNPsMat = NormalizedSNPsMat(I,:); NumSNPs = length(I); % changed number of SNPs (might be reduced in the intersection)
median_intensity_mat(:,2) = median(NormalizedSNPsMat)';
% Get Copy number : divide by normals average and multiply by 2
X_Inds = find(SNPChipAnnotStruct.chr_vec == 23); % Seperate treatment for X chromosome

[X_ids I_x J_x] = intersect(snp_ids, SNPChipAnnotStruct.snp_ids(X_Inds));
Y_Inds = find(SNPChipAnnotStruct.chr_vec == 24); % Seperate treatment for Y chromosome
[Y_ids I_y J_y] = intersect(snp_ids, SNPChipAnnotStruct.snp_ids(Y_Inds));
I_autosomal = setdiff(1:NumSNPs, union(I_x,I_y)); MaleInds = setdiff(1:NumSamples, FemaleInds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(std_iterate_flag)
    divide_cell = divide_by_normals_iterate(NormalizedSNPsMat(I_autosomal, :), NormalsInds);
    for i = 1:length(NormalizedSNPsMat, 2)
        NormalizedSNPsMat(I_autosomal,i) = 2*NormalizedSNPsMat(I_autosomal,i) ./ ...
            median(NormalizedSNPsMat(I_autosomal,divide_cell{i}),2);
        new_female_ind = intersect(FemaleInds, divide_cell{i});
        if(~length(new_female_ind))
            new_female_ind = FemaleInds;
        end
        NormalizedSNPsMat(I_x,i) = 2*NormalizedSNPsMat(I_x,i) ./ ...
            median(NormalizedSNPsMat(I_x,new_female_ind),2);

        new_male_ind = intersect(MaleInds, divide_cell{i});
        if(~length(new_male_ind))
            new_male_ind = MaleInds;
        end
        NormalizedSNPsMat(I_y,i) = NormalizedSNPsMat(I_y,i) ./ ...
            median(NormalizedSNPsMat(I_y,intersect(new_male_ind, NormalsInds)),2);
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oldNormalizedSNPsMat = NormalizedSNPsMat;
    NormalizedSNPsMat(I_autosomal,:) = 2*NormalizedSNPsMat(I_autosomal,:) ./ ...
        repmat(median(NormalizedSNPsMat(I_autosomal,NormalsInds),2), 1, NumSamples);
    NormalizedSNPsMat(I_x,:) = 2*NormalizedSNPsMat(I_x,:) ./ ...
        repmat(median(NormalizedSNPsMat(I_x,intersect(FemaleInds, NormalsInds)),2), 1, NumSamples);
    NormalizedSNPsMat(I_y,:) = NormalizedSNPsMat(I_y,:) ./ ...
        repmat(median(NormalizedSNPsMat(I_y,intersect(MaleInds, NormalsInds)),2), 1, NumSamples);
end

% calc % outliers in each array
array_outlier_vec = calc_array_outlier(NormalizedSNPsMat, std_outlier);
outlier_ind = find(array_outlier_vec>outlier_thresh);
if(remove_outlier_from_normals)
    NormalsInds = setdiff(NormalsInds, outlier_ind);
    if(length(NormalsInds))
        NormalizedSNPsMat(I_autosomal,:) = 2*oldNormalizedSNPsMat(I_autosomal,:) ./ ...
            repmat(median(oldNormalizedSNPsMat(I_autosomal,NormalsInds),2), 1, NumSamples);
        if(length(intersect(FemaleInds, NormalsInds)))
            NormalizedSNPsMat(I_x,:) = 2*oldNormalizedSNPsMat(I_x,:) ./ ...
                repmat(median(oldNormalizedSNPsMat(I_x,intersect(FemaleInds, NormalsInds)),2), 1, NumSamples);
        end
        if(length(intersect(MaleInds, NormalsInds)))
            NormalizedSNPsMat(I_y,:) = oldNormalizedSNPsMat(I_y,:) ./ ...
                repmat(median(oldNormalizedSNPsMat(I_y,intersect(MaleInds, NormalsInds)),2), 1, NumSamples);
        end
    end
end

clear oldNormalizedSNPsMat;

% Save in input format for the HMM
load(fullfile(user_dir, ['AlleleRatioMat_' lower(chip) '.mat']), 'allele_ratio_mat');

for i=1:NumSamples
    sample_str = [SampleNames{i} '_' lower(chip) '.mat'];
    copy_num_vec = NormalizedSNPsMat(:,i);
    allele_ratio_vec = allele_ratio_mat(:,i);
    save(fullfile(user_dir, sample_str), 'snp_ids', 'copy_num_vec', 'allele_ratio_vec');
end

% Addition: Save also one big matrix with everything. Save A and B intensities
[NormalizedSNPsCopyMatA, NormalizedSNPsCopyMatB] = RatioToCopyMats(NormalizedSNPsMat, allele_ratio_mat);

save(fullfile(user_dir, ['AllSamplesMat_' lower(chip) '.mat']), 'snp_ids', 'SampleNames', 'chip', ...
    'NormalizedSNPsCopyMatA', 'NormalizedSNPsCopyMatB', 'array_outlier_vec');

% create save normalizarion summary file
create_and_save_norm_summary_file(SampleNames, median_intensity_mat, user_dir, chip, array_outlier_vec);
