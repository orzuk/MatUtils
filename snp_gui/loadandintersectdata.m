% Written by Or Zuk 7/2007
%
% This function reads the SNPs data in the hmm input format from a .mat file and
% performs intersection with the SNPs which appear on the chip file and with
% the SNPs appearing in the HAPMAP database
%
% The inputs are:
%
% user_dir - directory with data
% sample_name - name of sample to read
% LDStruct - Structure containing Linkeage Disequilibrium information
% SNPChipAnnotStruct - Structure containing information on current chip
% HMMParamsStruct - Structure containing HMM parapeters
% HMMOutStructIn - previous output of hmm program
%
% The outputs are:
%
% ALL_MAT - Gives information on each SNP
% CHROM_MATS - The same (Gives information on each SNP) but divided to chromosomes
% HMMOutStruct - New output of hmm program
function [ALL_MAT CHROM_MATS HMMOutStruct] = LoadAndIntersectData(user_dir, sample_name, LDStruct, ...
    SNPChipAnnotStruct, HMMParamsStruct, HMMOutStructIn)

AssignAllGlobalConstants();

HMMOutStruct = HMMOutStructIn;
chip_type = SNPChipAnnotStruct.chip;
load([user_dir '\' sample_name '_' chip_type '.mat']);

% In the new version, we need to do some conversions of the data ...
% Work for now only on the current SNP chip
%% snp_id_str = ['snp_id_' lower(chip_type)];
%% eval(['snp_id_chip = ' snp_id_str ';']); % Just transfer to a convenient name ..
snp_id_chip = snp_ids;
[SnpsNames I J] = intersect(SNPChipAnnotStruct.snp_ids, snp_id_chip); % No relation to chromosome here
SNPChipAnnotStruct.rs_ids2 = SNPChipAnnotStruct.rs_ids(I); SNPChipAnnotStruct.snp_ids2 = SNPChipAnnotStruct.snp_ids(I); % Keep also the snp ids. Do not destroy the original vec.
HMM_locations = SNPChipAnnotStruct.chr_loc_vec(I); % Currently we don't know the locations
HMM_SNP_IDs = SNPChipAnnotStruct.snp_ids(I); % Save also the SNP IDs - why not ?
HMM_samples = {}; HMM_samples{1} = sample_name; % Pick one of the samples
HMM_ref_labels = 1;
HMM_chromosome_arm = SNPChipAnnotStruct.chr_vec(I); % Problem ! here we've got only the chromosome and not the arm!
allele_ratio_vec = min(allele_ratio_vec, 9999999999);

HMM_data_B = copy_num_vec ./ (allele_ratio_vec + 1); % We assume sample_ratio is A/B
HMM_data_A = HMM_data_B .* allele_ratio_vec;
HMM_data_B = HMM_data_B(J); HMM_data_A = HMM_data_A(J);


if(HMMParamsStruct.use_affy_genotypes) % Here we know the Affymetrix genotypes (doesn't make sense unless we debug stuff)
    HMM_data_geno = genotype_vec;
    HMM_data_geno = HMM_data_geno(J); % Sort
end

epsilon = 0.000000001; 
do_log = 0; % Flag saying if to perform log transform
place_flag = 1; % give a different distribution for every place. These are based on chromosome 3
% which is divided to 10 patients with LOH and 10 patients with NLOH

num_genes = length(snp_id_chip); num_samples = length(HMM_samples);

% First do fold change if neccessary. We will keep also the healthy ones !!!!!!
do_fold_change=0;
if(do_fold_change)
    HMM_ref_mean = mean(HMM_data_B(:,HMM_ref_labels==0), 2);  % Must be first !!!
    num_samples = length(HMM_samples); % Correct the number of samples
    % Now make the fold-change !!! and take the log !!
    HMM_data_B = HMM_data_B ./ repmat(HMM_ref_mean, 1, num_samples);
    HMM_data_B = log(HMM_data_B);
    % Make a plot of the means of patients log-foLDStruct.LD-change at each chromosome
end

% Get rid of this stupid cell format, and get data as numbers
HMM_chromosome = HMM_chromosome_arm; % Here we already have the chromosomes only
%     HMM_chromosome = (char(HMM_chromosome_arm)); % Here we already have the chromosomes only
%     XXX = find(HMM_chromosome == 'X');
%     HMM_chromosome(XXX,1) = '2'; % Encode the X chromosome
%     HMM_chromosome(XXX,2) = '3'; % Encode the X chromosome
%     III = intersect( find(HMM_chromosome(:,1) == ' '), find(HMM_chromosome(:,2) == ' ') );
%     HMM_chromosome(III,1) = '0'; % Encode the chromosomes not known at all
%     HMM_chromosome = str2num(HMM_chromosome);

CHROM_MATS.loc = {}; CHROM_MATS.data_A = {}; CHROM_MATS.data_B = {}; % Now we set the data also in cell file !!
ALL_MAT.data = []; ALL_MAT.data_A = []; ALL_MAT.data_B = [];
ALL_MAT.loc = []; ALL_MAT.snp_ids = []; ALL_MAT.strand = [];
if(HMMParamsStruct.use_affy_genotypes)
    ALL_MAT.data_geno = [];
end

for i=HMMParamsStruct.ChromosomesToRun  % loop over all chromosomes

    indexes_vec = find(HMM_chromosome == i);     % find which SNPs are on current chromosome
    [JointSNPs III JJJ] = intersect(LDStruct.LD{i}.RS_IDs, SNPChipAnnotStruct.rs_ids2(indexes_vec)); % need to do further intersection since some SNPs are ? have the same SNPChipAnnotStruct.rs_ids2 ?
    indexes_vec = indexes_vec(JJJ);
    MissInds = setdiff( [1:length(LDStruct.LD{i}.RS_IDs)], III );
    HMMOutStruct.SNPsIDs{i}(III) = SNPChipAnnotStruct.snp_ids2(indexes_vec);
    for j=MissInds  % We add these guys which do not have a SNP i.d.
        HMMOutStruct.SNPsIDs{i}{j} = ' ';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Debug - serch for empty SNPs
    for j=1:length(SNPChipAnnotStruct.snp_ids2)
        if(isempty(SNPChipAnnotStruct.snp_ids2{j}))
            empty_j = j
        end
    end

    if(~isempty(LDStruct.LD{i}.RS_IDs))
        % Take the sub-matrix of SNPs which appear on the current chromsome
        CHROM_MATS.data_A{i} = zeros(length(LDStruct.LD{i}.RS_IDs),1);
        CHROM_MATS.data_B{i} = zeros(length(LDStruct.LD{i}.RS_IDs),1);
        CHROM_MATS.loc{i} = zeros(length(LDStruct.LD{i}.RS_IDs),1);
        CHROM_MATS.data_A{i}(III) = HMM_data_A(indexes_vec,:);
        CHROM_MATS.data_B{i}(III) = HMM_data_B(indexes_vec,:);

        CHROM_MATS.loc{i}(III) = HMM_locations(indexes_vec);
        CHROM_MATS.snp_ids{i}{length(LDStruct.LD{i}.RS_IDs)} = '';
        CHROM_MATS.snp_ids{i}(III) = HMM_SNP_IDs(indexes_vec);
        if(HMMParamsStruct.use_affy_genotypes)
            CHROM_MATS.data_geno{i} = zeros(length(LDStruct.LD{i}.RS_IDs),1); % Add also Affy genotypes for 'sanity check'
            CHROM_MATS.data_geno{i}(III) = HMM_data_geno(indexes_vec,:);
        end

        for j=MissInds % These are not in the intersection
            CHROM_MATS.data_A{i}(j) = ...
                (CHROM_MATS.data_A{i}(max(j-1,1)) + CHROM_MATS.data_A{i}(min(j+1,length(LDStruct.LD{i}.RS_IDs))))/2;
            CHROM_MATS.data_B{i}(j) = ...
                (CHROM_MATS.data_B{i}(max(j-1,1)) + CHROM_MATS.data_B{i}(min(j+1,length(LDStruct.LD{i}.RS_IDs))))/2;
            if(HMMParamsStruct.use_affy_genotypes)
                CHROM_MATS.data_geno{i}(j) = 4; % unknown genotype
            end
            CHROM_MATS.loc{i}(j) = ...
                (CHROM_MATS.loc{i}(max(j-1,1)) + CHROM_MATS.loc{i}(min(j+1,length(LDStruct.LD{i}.RS_IDs))))/2;
            CHROM_MATS.snp_ids{i}{j} = ''; % SNP unknown
        end

        % Calc the 'singleton' LDStruct.LD transition matrices. This means that we take into account only frequencies,
        % and NOT the pairwise LDStruct.LD correlations
        HT_SIGNS = zeros(1,length(I));
        HT_SIGNS(strmatch('+', SNPChipAnnotStruct.strand(I))) = 1; % '+' SNPChipAnnotStruct.strand is 1
        CHROM_MATS.strand{i} = zeros(1,length(LDStruct.LD{i}.FreqVec));
        CHROM_MATS.strand{i}(III) = HT_SIGNS(indexes_vec);
        %%% LDStruct.LD_corrected_FreqVec = 2 .* LDStruct.LD{i}.FreqVec .* CHR_HT_SIGNS' - LDStruct.LD{4}.FreqVec - CHR_HT_SIGNS' + 1;

        % Change data according to strand (why do it here?)
        tmp_Y = CHROM_MATS.data_A{i};
        CHROM_MATS.data_A{i} = CHROM_MATS.data_A{i} .* CHROM_MATS.strand{i}' + CHROM_MATS.data_B{i} .* (1-CHROM_MATS.strand{i}');
        CHROM_MATS.data_B{i} = tmp_Y .* (1-CHROM_MATS.strand{i}') + CHROM_MATS.data_B{i} .* CHROM_MATS.strand{i}'; % Switch Y and Y2 where SNPChipAnnotStruct.strand is negative
        if(HMMParamsStruct.use_affy_genotypes)
            CHROM_MATS.data_geno{i} =  CHROM_MATS.data_geno{i} + ...
                mod(CHROM_MATS.data_geno{i},2) .* (1-CHROM_MATS.strand{i}') .* ...
                2 .* (2-CHROM_MATS.data_geno{i}); % Switch 1 AA with 3 BB
        end

        % Sort data according to chromosomal locations
        [ CHROM_MATS.loc{i} perm ] = sort(CHROM_MATS.loc{i});
        CHROM_MATS.data_A{i} = CHROM_MATS.data_A{i}(perm);
        CHROM_MATS.data_B{i} = CHROM_MATS.data_B{i}(perm);
        CHROM_MATS.strand{i} = CHROM_MATS.strand{i}(perm);
        CHROM_MATS.snp_ids{i} = CHROM_MATS.snp_ids{i}(perm); % Sort also the SNP IDs
        HMMOutStruct.SNPsIDs{i} = HMMOutStruct.SNPsIDs{i}(perm);

        % Concatenate all data from all chromosomes
        ALL_MAT.data = [ALL_MAT.data CHROM_MATS.data_A{i}'+CHROM_MATS.data_B{i}'];
        ALL_MAT.data_A = [ALL_MAT.data_A CHROM_MATS.data_A{i}'];
        ALL_MAT.data_B = [ALL_MAT.data_B CHROM_MATS.data_B{i}'];
        ALL_MAT.strand = [ALL_MAT.strand CHROM_MATS.strand{i}]; % No transpose here ...
        ALL_MAT.loc = [ALL_MAT.loc CHROM_MATS.loc{i}']; % This should be changed to account for different chromosomes
        ALL_MAT.snp_ids = [ALL_MAT.snp_ids CHROM_MATS.snp_ids{i}]; % This should be changed to account for different chromosomes
        ALL_MAT.chr = i; % Save the chromosomes

        if(HMMParamsStruct.use_affy_genotypes)
            CHROM_MATS.data_geno{i} = CHROM_MATS.data_geno{i}(perm);
            ALL_MAT.data_geno = [ALL_MAT.data_geno CHROM_MATS.data_geno{i}'];
        end

    end
end

sprintf('Finished Collecting Data')