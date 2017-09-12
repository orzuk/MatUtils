% Note : We return here SINGLE copy and INT genotype - To save Memory!!! 
function [data_A data_B genotypes_AB snp_ids strand chr] = UpdateOnePopulation(hapmap_population, user_dir, chip_type, ...
    SampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct, hapmap_dir, hapmap_version)

ALL_MAT = {};
AssignAllGlobalConstants();
HMMParamsStruct.use_affy_genotypes = 0; % Should be eliminated - input as a button ...
chip_type = lower(SNPChipAnnotStruct.chip);
num_samples = length(SampleNames);
chroms = 1:24;

sample_by_sample_flag = 0; 

if(sample_by_sample_flag)
    for(cur_sample = 1:num_samples) % Outer loop on samples .... create one file for all samples
        sample_name = SampleNames{cur_sample}
        HMMOutStruct = {};
        for i=chroms
            HMMOutStruct.SNPsIDs{i} = {};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Just make all the indexology right here
        [TMP_ALL_MAT CHROM_MATS HMMOutStruct] = LoadAndIntersectData(user_dir, sample_name, LDStruct, SNPChipAnnotStruct, ...
            HMMParamsStruct, HMMOutStruct);
        num_snps = length(TMP_ALL_MAT.loc);

        if(cur_sample==1) % On the first sample prepare the matrices
            snp_ids = TMP_ALL_MAT.snp_ids;
            ALL_MAT.loc = TMP_ALL_MAT.loc;
            strand = TMP_ALL_MAT.strand;
            chr = TMP_ALL_MAT.chr;

            data_A = zeros(num_samples, num_snps, 'single');
            data_B = zeros(num_samples, num_snps, 'single');
        end

        data_A(cur_sample,:) = single(TMP_ALL_MAT.data_A);
        data_B(cur_sample,:) = single(TMP_ALL_MAT.data_B);
    end % loop for ???
else  % Take all samples from one file ..
    load(fullfile(user_dir, ['AllSamplesMat_', chip_type,  '.mat']));
    data_A = single(NormalizedSNPsCopyMatA'); data_B = single(NormalizedSNPsCopyMatB');
    clear NormalizedSNPsCopyMatA NormalizedSNPsCopyMatB;
    ChipStruct = load(fullfile('../database', [chip_type '_annot_data_' genome_assembly '.mat']), 'snp_ids', 'strand', 'chr_vec');
    [S I J] = intersect(snp_ids, ChipStruct.snp_ids);
    strand = cell(length(snp_ids),1);
    strand(:) = {''};
    strand(I) = ChipStruct.strand(J);
    chr = zeros(length(snp_ids),1);
    chr(I) = ChipStruct.chr_vec(J);
end


clear CHROM_MATS HMMOutStruct;

SampleInds = GetHapMapSampleIndex(SampleNames, hapmap_population); SampleBits = mod(SampleInds-1,30)+1; SampleWords = ceil(SampleInds./30);

% load the HAPMAP SNPs which appear on the chip (e.g. hind, xba or both)
R=load( fullfile(hapmap_dir, pop_str_vec{hapmap_population}, hapmap_version, ...
    [chip_type, '_genotypes_chr_' pop_str_vec{hapmap_population} '_' hapmap_version '_nr_fwd.mat']) );

genotypes_AB = zeros(num_samples, length(snp_ids), 'uint8');
for chrom=chroms
    chrom_is = chrom
    [snp_ids_intersect I J] = intersect(snp_ids, R.SnpsIDs{chrom}); % Get the common SNPs from disp

    CurSNPsDataChromMat = R.SnpsData{chrom}(J,:);
    CurSNPsBadCallChrom = find(bitget(R.SnpsBadCalls{chrom}(J,1),1));

    for cur_sample = 1:num_samples
        sample_is = cur_sample

        HAPMAP_SNPs = 2*bitget(CurSNPsDataChromMat(:,SampleWords(cur_sample)), SampleBits(cur_sample)) + ...
            bitget(CurSNPsDataChromMat(:,SampleWords(cur_sample)+3), SampleBits(cur_sample)); % Take the first patient from HAPMAP. Bad Calls are considered as '-1'
        HAPMAP_SNPs(CurSNPsBadCallChrom) = -1; % Bad calls NN SNPs are set to minus one

        HAPMAP_SNPs(HAPMAP_SNPs == 2) = 1; % unite AB and BA

        % Transfer from (0,1,3,-1) to (1,2,3,4)
        HAPMAP_SNPs(HAPMAP_SNPs == 3) = BB;
        HAPMAP_SNPs(HAPMAP_SNPs == 1) = AB;
        HAPMAP_SNPs(HAPMAP_SNPs == 0) = AA;
        HAPMAP_SNPs(HAPMAP_SNPs == -1) = NoCall;

        genotypes_AB(cur_sample,I) = uint8(HAPMAP_SNPs);
    end
end

