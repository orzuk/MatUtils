% Written by Or Zuk 5/2007
%
% This function reads the genotypes derived from the HAPMAP data, and
% compares it to the output of the HMM program, which is derived from the
% SNP-chip data that we analyzed. Clearly, one cannot expect to get 100%
% accuracy. A more reasonable goal is to get comparable, or better results,
% then those obtained by Affymetrix genotyping. To measure this, we also
% compare here the affymetrix genotypes with the HAPMAP genotypes, and
% obtain their accuracy.
% 
% Input:
% SampleNames - names of samples to compare
% chip_type - type of chip
% user_dir - working directory
% comp_hap_affy - flag saying if to compare hapmap and affymetrix genotypes
% comp_hap_hmm - flag saying if to compare hapmap and hmm genotypes
% comp_affy_hmm - flag saying if to compare affymetrix and hmm genotypes
% use_strand_vec - vector of flags saying whether or not to reverse strand in comparison.
% hapmap_dir - directory where the hapmap stuff lies
% hapmap_population - which population do we deal with
% hapmap_version - version of hapmap data
% 
% Output: 
% ErrorStruct - containing all the error parameters
%
function ErrorStruct = HapMapCompareAffyAndHMMGenotypes(num_samples, chip_type, user_dir, chroms, ...
    comp_hap_affy, comp_hap_hmm, comp_affy_hmm, use_strand_vec, hapmap_dir, hapmap_population, hapmap_version, RLMM_population)

AssignAllGlobalConstants();
num_chroms = length(chroms);
%chroms = 1:22; % Currently do all the chromosomes
do_couples=1; joint_flag = 1; % Flags for getting marginals for HMM

   
display_dir = ['display_' pop_str_vec{RLMM_population}];
load(fullfile(user_dir, display_dir, ['all_samples_hmm_out_' chip_type '.mat']),...
    'sample_names', 'data_snp_ids', 'genotype_mat');

if(num_samples == -1)
    nsamples = length(sample_names);
else
    nsamples = min(num_samples, length(sample_names));
end

comp_hap_local = 0; % ~isempty(data_snp_ids); currently not working

if(comp_hap_local)
    ErrorStruct.hap_local_chr_common_snps_vec = zeros(1,num_chroms); % How many SNPs are on the chip and HAPMAP in each chrom
    ErrorStruct.hap_local_chr_error_snps_vec = zeros(nsamples,num_chroms); % For how many SNPs did we get a different genotype
    ErrorStruct.hap_local_uneq_inds_mat = zeros(nsamples, length(data_snp_ids));
    % Special temp
    ErrorStruct.hap_nocalls_mat = zeros(nsamples, length(data_snp_ids));
end

if(comp_hap_hmm)
    ErrorStruct.hap_hmm_chr_common_snps_vec = zeros(1,num_chroms); % How many SNPs are on the chip and HAPMAP in each chrom
    ErrorStruct.hap_hmm_chr_error_snps_vec = zeros(nsamples,num_chroms); % For how many SNPs did we get a different genotype
end
if(comp_hap_affy)
    ErrorStruct.hap_affy_chr_common_snps_vec = zeros(1,num_chroms); % How many SNPs are on the chip and HAPMAP in each chrom for Affymetrix
    ErrorStruct.hap_affy_chr_error_snps_vec = zeros(nsamples,num_chroms); % For how many SNPs did we get a different genotype for Affymetrix
    ErrorStruct.affy_no_calls = zeros(nsamples,num_chroms);
end
if(comp_affy_hmm)
    ErrorStruct.affy_hmm_chr_common_snps_vec = zeros(1,num_chroms); % How many SNPs are on the chip and HMM in each chrom for Affymetrix
    ErrorStruct.affy_hmm_chr_error_snps_vec = zeros(nsamples,num_chroms); % For how many SNPs did we get a different genotype for Affymetrix
    ErrorStruct.affy_no_calls = zeros(nsamples,num_chroms);
end

% load the HAPMAP SNPs which appear on the chip (e.g. hind, xba or both) ...
if(comp_hap_affy || comp_hap_hmm)
    R=load(fullfile(hapmap_dir, pop_str_vec{hapmap_population}, hapmap_version, ...
        [chip_type, '_genotypes_chr_' pop_str_vec{hapmap_population} '_' hapmap_version '_nr_fwd.mat']));    
end

SampleInds = GetHapMapSampleIndex(sample_names, hapmap_population); SampleBits = mod(SampleInds-1,30)+1; SampleWords = ceil(SampleInds./30);

% Load Chip Annotations - We need the strand from here !!!
SNPChipAnnotStruct = load(['..\database\' chip_type '_annot_data_' genome_assembly '.mat'], 'snp_ids', 'strand');


% Here we need to reverse the relevant SNPs according to strand
StrandSigns = zeros(1,length(SNPChipAnnotStruct.strand));
StrandSigns(strmatch('-', SNPChipAnnotStruct.strand)) = 1;

% New: Make the intersection before to gain speed
if(comp_hap_local)
    for chrom=chroms
        [snp_ids_intersect{chrom} III{chrom} JJJ{chrom}] = intersect(data_snp_ids, R.SnpsIDs{chrom}); % Get the common SNPs from disp
        ErrorStruct.hap_local_snp_ids(III{chrom}) = snp_ids_intersect{chrom};
    end
end

% Loop over all the samples
for cur_sample = 1:nsamples
    do_sample = cur_sample

    % Load the results of the Affymetrix genotypes
    if(comp_hap_affy || comp_affy_hmm)
        % Should we load also the allele ratios? Yes if we want to plot stuff!
        Aff = load(fullfile(user_dir, [sample_names{cur_sample} '_' chip_type '.mat']));  
%        eval(['Aff = load(fullfile(user_dir, [sample_names{cur_sample} ''_'' chip_type ''.mat'']));']);  % Should we load also the allele ratios? Yes if we want to plot stuff!
    end

    % Load the results of the HMM - from the display !! 
    if(comp_hap_hmm || comp_affy_hmm)
        load(fullfile(user_dir, display_dir , [sample_names{cur_sample} '_' chip_type '_disp']));
    end


    % Set the relevant string and strand - NOTE: We do it for the HAPMAP data and NOT for the affy and/or HMM genotypes !!!
    if(comp_hap_affy || comp_affy_hmm)
        [snp_ids_intersect I J] = intersect(SNPChipAnnotStruct.snp_ids, Aff.snp_ids);
    end
    %%     eval(['Aff.genotype_vec_' chip_type '(J(find(StrandSigns(I)))) = 4-Aff.genotype_vec_' chip_type '(J(find(StrandSigns(I))));']);
    %%     eval(['Aff.genotype_vec_' chip_type '(find(Aff.genotype_vec_' chip_type ' == 0)) = 4;']);

    % Go chrom-by-chrom as to not take too much memory
    ErrorStruct.hap_local_ErrorMat{cur_sample} = zeros(4);
    ErrorStruct.hap_affy_ErrorMat{cur_sample} = zeros(4);
    ErrorStruct.hap_hmm_ErrorMat{cur_sample} = zeros(4);
    ErrorStruct.affy_hmm_ErrorMat{cur_sample} = zeros(4);

    for chrom=chroms
        do_chrom = chrom

        if(comp_hap_local)
%%            [snp_ids_intersect I J] = intersect(data_snp_ids, R.SnpsIDs{chrom}); % Get the common SNPs from disp
            ErrorStruct.hap_local_chr_common_snps_vec(chrom) = length(III{chrom});
            HAPMAP_SNPs = 2*bitget(R.SnpsData{chrom}(JJJ{chrom},SampleWords(cur_sample)), SampleBits(cur_sample)) + ...
                bitget(R.SnpsData{chrom}(JJJ{chrom},SampleWords(cur_sample)+3), SampleBits(cur_sample)); % Take the first patient from HAPMAP. Bad Calls are considered as '-1'
            if(use_strand_vec(4))             % Correct for reverse strand
                [strand_snp_ids_intersect strand_I strand_J] = intersect(SNPChipAnnotStruct.snp_ids,  R.SnpsIDs{chrom}(JJJ{chrom})); % Get the common SNPs
                HAPMAP_SNPs(strand_J(find(StrandSigns(strand_I)))) = 3-HAPMAP_SNPs(strand_J(find(StrandSigns(strand_I))));
            end
            HAPMAP_SNPs(HAPMAP_SNPs ==2) = 1;

            HAPMAP_SNPs(find(bitget(R.SnpsBadCalls{chrom}(JJJ{chrom},1), 1))) = -1; % Bad calls NN SNPs are set to minus one
            
            LOCAL_SNPs = genotype_mat(III{chrom},cur_sample);
            LOCAL_SNPs(LOCAL_SNPs == AA) = 0; % Transfer (1,2,3) to (0,1,3)           
            LOCAL_SNPs(LOCAL_SNPs == AB) = 1; % AB and BA are the same ..
            LOCAL_SNPs(LOCAL_SNPs == BB) = 3; % Transfer (1,2,3) to (0,1,3)        
            
            ErrorStruct.hap_local_chr_error_snps_vec(cur_sample,chrom) = sum(LOCAL_SNPs ~= HAPMAP_SNPs); % Get the number of errors
            ErrorStruct.hap_local_uneq_inds_mat(cur_sample, III{chrom}) = (LOCAL_SNPs ~= HAPMAP_SNPs); % Set the matrix
            ErrorStruct.hap_nocalls_mat(cur_sample, III{chrom}) = (HAPMAP_SNPs == -1); % Set the bad calls 
            
            % Generate also a matrix to see where the mistakes are .. 
            ErrorStruct.hap_local_ErrorMat{cur_sample}(1,2) = ErrorStruct.hap_local_ErrorMat{cur_sample}(1,2) + length(intersect( find(HAPMAP_SNPs == -1), find(LOCAL_SNPs == 0) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(1,3) = ErrorStruct.hap_local_ErrorMat{cur_sample}(1,3) + length(intersect( find(HAPMAP_SNPs == -1), find(LOCAL_SNPs == 1) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(1,4) = ErrorStruct.hap_local_ErrorMat{cur_sample}(1,4) + length(intersect( find(HAPMAP_SNPs == -1), find(LOCAL_SNPs == 3) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(2,3) = ErrorStruct.hap_local_ErrorMat{cur_sample}(2,3) + length(intersect( find(HAPMAP_SNPs == 0), find(LOCAL_SNPs == 1) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(2,4) = ErrorStruct.hap_local_ErrorMat{cur_sample}(2,4) + length(intersect( find(HAPMAP_SNPs == 0), find(LOCAL_SNPs == 3) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(3,2) = ErrorStruct.hap_local_ErrorMat{cur_sample}(3,2) + length(intersect( find(HAPMAP_SNPs == 1), find(LOCAL_SNPs == 0) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(3,4) = ErrorStruct.hap_local_ErrorMat{cur_sample}(3,4) + length(intersect( find(HAPMAP_SNPs == 1), find(LOCAL_SNPs == 3) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(4,2) = ErrorStruct.hap_local_ErrorMat{cur_sample}(4,2) + length(intersect( find(HAPMAP_SNPs == 3), find(LOCAL_SNPs == 0) ));
            ErrorStruct.hap_local_ErrorMat{cur_sample}(4,3) = ErrorStruct.hap_local_ErrorMat{cur_sample}(4,3) + length(intersect( find(HAPMAP_SNPs == 3), find(LOCAL_SNPs == 1) ));
        
        
        end
        
        
        if(comp_hap_hmm)
            %[snp_ids_intersect I J] = intersect(DispStruct.Chrom{chrom}.SNPsIDs, R.SnpsIDs{chrom}); % Get the common SNPs from disp
            [snp_ids_intersect I J] = intersect(data_snp_ids, R.SnpsIDs{chrom}); % Get the common SNPs from disp
            
            ErrorStruct.hap_hmm_chr_common_snps_vec(chrom) = length(I);


            HAPMAP_SNPs = 2*bitget(R.SnpsData{chrom}(J,SampleWords(cur_sample)), SampleBits(cur_sample)) + ...
                bitget(R.SnpsData{chrom}(J,SampleWords(cur_sample)+3), SampleBits(cur_sample)); % Take the first patient from HAPMAP. Bad Calls are considered as '-1'
            HAPMAP_SNPs(find(bitget(R.SnpsBadCalls{chrom}(J,1), 1))) = -1; % Bad calls NN SNPs are set to minus one

            % New: We get the input from the display instead of the hmm_out
            %HMMOut_SNPs = DispStruct.Chrom{chrom}.Genotypes(I);  % Take the patient from the HMM Output display file
            HMMOut_SNPs = genotype_mat(I, cur_sample);  % Take the patient from the HMM Output display file
            
            %%            HMMOut_SNPs = VitStruct.joint_genotype(I);  % Take the patient from the HMM Output ..
            HMMOut_SNPs(HMMOut_SNPs == 2) = 1; % AB and BA are the same ..
            ErrorStruct.hap_hmm_chr_error_snps_vec(cur_sample,chrom) = sum(HMMOut_SNPs ~= HAPMAP_SNPs); % Get the number of errors
            %%      figure; plot(HMMOut_SNPs+0.1*randn(length(HMMOut_SNPs),1), HAPMAP_SNPs+0.1*randn(length(HAPMAP_SNPs),1), '.'); xlabel('HMM'); ylabel('HapMap');

        
            % Generate also a matrix to see where the mistakes are .. 
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(1,2) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(1,2) + length(intersect( find(HAPMAP_SNPs == -1), find(HMMOut_SNPs == 0) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(1,3) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(1,3) + length(intersect( find(HAPMAP_SNPs == -1), find(HMMOut_SNPs == 1) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(1,4) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(1,4) + length(intersect( find(HAPMAP_SNPs == -1), find(HMMOut_SNPs == 3) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(2,3) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(2,3) + length(intersect( find(HAPMAP_SNPs == 0), find(HMMOut_SNPs == 1) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(2,4) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(2,4) + length(intersect( find(HAPMAP_SNPs == 0), find(HMMOut_SNPs == 3) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(3,2) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(3,2) + length(intersect( find(HAPMAP_SNPs == 1), find(HMMOut_SNPs == 0) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(3,4) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(3,4) + length(intersect( find(HAPMAP_SNPs == 1), find(HMMOut_SNPs == 3) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(4,2) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(4,2) + length(intersect( find(HAPMAP_SNPs == 3), find(HMMOut_SNPs == 0) ));
            ErrorStruct.hap_hmm_ErrorMat{cur_sample}(4,3) = ErrorStruct.hap_hmm_ErrorMat{cur_sample}(4,3) + length(intersect( find(HAPMAP_SNPs == 3), find(HMMOut_SNPs == 1) ));
        
        end

        if(comp_hap_affy)
            [hap_affy_snp_ids_intersect hap_affy_I hap_affy_J] = intersect(Aff.snp_ids, R.SnpsIDs{chrom}); % Get the common SNPs
            ErrorStruct.hap_affy_chr_common_snps_vec(chrom) = length(hap_affy_I);

            HAPMAP_SNPs = 2*bitget(R.SnpsData{chrom}(hap_affy_J,SampleWords(cur_sample)), SampleBits(cur_sample)) + ...
                bitget(R.SnpsData{chrom}(hap_affy_J,SampleWords(cur_sample)+3), SampleBits(cur_sample)); % Take the first patient from HAPMAP
            if(use_strand_vec(2))             % Correct for reverse strand
                [strand_snp_ids_intersect strand_I strand_J] = intersect(SNPChipAnnotStruct.snp_ids,  R.SnpsIDs{chrom}(hap_affy_J)); % Get the common SNPs
                HAPMAP_SNPs(strand_J(find(StrandSigns(strand_I)))) = 3-HAPMAP_SNPs(strand_J(find(StrandSigns(strand_I))));
            end
            HAPMAP_SNPs(HAPMAP_SNPs ==2) = 1;
            HAPMAP_SNPs(find(bitget(R.SnpsBadCalls{chrom}(hap_affy_J,1), 1))) = -1; % Bad calls NN SNPs are set to minus one

            AFFY_SNPs = Aff.genotype_vec(hap_affy_I)-1; % Get the SNPs that were obtained by affymetrix software from their chips ...
            AFFY_SNPs(AFFY_SNPs == 3) = -1; % The No-Calls
            AFFY_SNPs(AFFY_SNPs == 2) = 3; % The BB's
            ErrorStruct.hap_affy_chr_error_snps_vec(cur_sample,chrom) = sum(AFFY_SNPs ~= HAPMAP_SNPs); % Get the number of errors
            ErrorStruct.affy_no_calls(cur_sample,chrom) = sum(AFFY_SNPs == -1);

            % Generate also a matrix to see where the mistakes are .. 
            
            %%            ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,1) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,1) + length(intersect( find(AFFY_SNPs == -1), find(HAPMAP_SNPs == -1) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,2) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,2) + length(intersect( find(HAPMAP_SNPs == -1), find(AFFY_SNPs == 0) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,3) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,3) + length(intersect( find(HAPMAP_SNPs == -1), find(AFFY_SNPs == 1) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,4) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(1,4) + length(intersect( find(HAPMAP_SNPs == -1), find(AFFY_SNPs == 3) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(2,1) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(2,1) + length(intersect( find(HAPMAP_SNPs == 0), find(AFFY_SNPs == -1) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(2,3) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(2,3) + length(intersect( find(HAPMAP_SNPs == 0), find(AFFY_SNPs == 1) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(2,4) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(2,4) + length(intersect( find(HAPMAP_SNPs == 0), find(AFFY_SNPs == 3) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(3,1) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(3,1) + length(intersect( find(HAPMAP_SNPs == 1), find(AFFY_SNPs == -1) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(3,2) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(3,2) + length(intersect( find(HAPMAP_SNPs == 1), find(AFFY_SNPs == 0) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(3,4) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(3,4) + length(intersect( find(HAPMAP_SNPs == 1), find(AFFY_SNPs == 3) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(4,1) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(4,1) + length(intersect( find(HAPMAP_SNPs == 3), find(AFFY_SNPs == -1) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(4,2) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(4,2) + length(intersect( find(HAPMAP_SNPs == 3), find(AFFY_SNPs == 0) ));
            ErrorStruct.hap_affy_ErrorMat{cur_sample}(4,3) = ErrorStruct.hap_affy_ErrorMat{cur_sample}(4,3) + length(intersect( find(HAPMAP_SNPs == 3), find(AFFY_SNPs == 1) ));
            %%        figure; plot(AFFY_SNPs+0.1*randn(length(AFFY_SNPs),1), HAPMAP_SNPs+0.1*randn(length(AFFY_SNPs),1), '.'); xlabel('Affy'); ylabel('HapMap');
            
        end

        if(comp_affy_hmm)
            [affy_hmm_snp_ids_intersect affy_hmm_I affy_hmm_J] = intersect(Aff.snp_ids, DispStruct.Chrom{chrom}.SNPsIDs); % Get the common SNPs

            ErrorStruct.affy_hmm_chr_common_snps_vec(chrom) = length(affy_hmm_I);
            AFFY_SNPs = Aff.genotype_vec(affy_hmm_I)-1; % Get the SNPs that were obtained by affymetrix software from their chips ...
            AFFY_SNPs(AFFY_SNPs == 3) = -1; % The No-Calls
            AFFY_SNPs(AFFY_SNPs == 2) = 3; % The BB's
            % New: We get the input from the display instead of the hmm_out
            HMMOut_SNPs = DispStruct.Chrom{chrom}.Genotypes(affy_hmm_J);  % Take the patient from the HMM Output display file
            HMMOut_SNPs(HMMOut_SNPs == 2) = 1; % AB and BA are the same ..

            % Correct for reverse strand - should we ? YES
            if(use_strand_vec(3))
                [snp_ids_intersect strand_I strand_J] = intersect(SNPChipAnnotStruct.snp_ids, DispStruct.Chrom{chrom}.SNPsIDs(affy_hmm_J));
                HMMOut_SNPs(strand_J(find(StrandSigns(strand_I)))) = 3-HMMOut_SNPs(strand_J(find(StrandSigns(strand_I))));
            end
            HMMOut_SNPs(HMMOut_SNPs == 2) = 1; % AB and BA are the same ..
            ErrorStruct.affy_hmm_chr_error_snps_vec(cur_sample,chrom) = sum(AFFY_SNPs ~= HMMOut_SNPs); % Get the number of errors
            ErrorStruct.affy_no_calls(cur_sample,chrom) = sum(AFFY_SNPs == -1);
            
            
            % Generate also a matrix to see where the mistakes are ..
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(1,2) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(1,2) + length(intersect( find(AFFY_SNPs == -1), find(HMMOut_SNPs == 0) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(1,3) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(1,3) + length(intersect( find(AFFY_SNPs == -1), find(HMMOut_SNPs == 1) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(1,4) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(1,4) + length(intersect( find(AFFY_SNPs == -1), find(HMMOut_SNPs == 3) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(2,3) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(2,3) + length(intersect( find(AFFY_SNPs == 0), find(HMMOut_SNPs == 1) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(2,4) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(2,4) + length(intersect( find(AFFY_SNPs == 0), find(HMMOut_SNPs == 3) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(3,2) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(3,2) + length(intersect( find(AFFY_SNPs == 1), find(HMMOut_SNPs == 0) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(3,4) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(3,4) + length(intersect( find(AFFY_SNPs == 1), find(HMMOut_SNPs == 3) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(4,2) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(4,2) + length(intersect( find(AFFY_SNPs == 3), find(HMMOut_SNPs == 0) ));
            ErrorStruct.affy_hmm_ErrorMat{cur_sample}(4,3) = ErrorStruct.affy_hmm_ErrorMat{cur_sample}(4,3) + length(intersect( find(AFFY_SNPs == 3), find(HMMOut_SNPs == 1) ));
            
        end

        if(chrom*cur_sample == 0) %% 1) %% plot or not? 
            if(comp_hap_affy)
                PlotAlleleRatios(Aff.copy_num_vec(hap_affy_I), Aff.allele_ratio_vec(hap_affy_I), max(0,HAPMAP_SNPs), ['HAPMAP genotypes chrom ' num2str(chrom) ' Sample ' sample_names{cur_sample}]);
                PlotAlleleRatios(Aff.copy_num_vec(hap_affy_I), Aff.allele_ratio_vec(hap_affy_I), Aff.genotype_vec(hap_affy_I), ['affy genotypes chrom ' num2str(chrom) ' Sample ' sample_names{cur_sample}]);
                PlotAlleleRatios(Aff.copy_num_vec(hap_affy_I), Aff.allele_ratio_vec(hap_affy_I), AFFY_SNPs == HAPMAP_SNPs, ['affy vs. hap genotypes  Sample ' sample_names{cur_sample}]);
            end

            if(comp_affy_hmm || comp_hap_affy)
                PlotAlleleRatios(Aff.copy_num_vec, Aff.allele_ratio_vec, Aff.genotype_vec, ['affy genotypes  Sample ' sample_names{cur_sample}]);
            end
           if(comp_hap_hmm)
               PlotAlleleRatios(Aff.copy_num_vec(affy_hmm_I), Aff.allele_ratio_vec(affy_hmm_I), HMMOut_SNPs == HAPMAP_SNPs, ['hmm genotypes vs. HapMap genotypes chrom ' num2str(chrom) ' Sample ' sample_names{cur_sample}]);
           end
            
            if(comp_affy_hmm)
                PlotAlleleRatios(Aff.copy_num_vec(affy_hmm_I), Aff.allele_ratio_vec(affy_hmm_I), max(0, AFFY_SNPs), ['affy genotypes chrom ' num2str(chrom) ' Sample ' sample_names{cur_sample}]); % Aff.genotype_vec(affy_hmm_I)
                PlotAlleleRatios(Aff.copy_num_vec(affy_hmm_I), Aff.allele_ratio_vec(affy_hmm_I), HMMOut_SNPs, ['hmm genotypes chrom ' num2str(chrom) ' Sample ' sample_names{cur_sample}]);
                PlotAlleleRatios(Aff.copy_num_vec(affy_hmm_I), Aff.allele_ratio_vec(affy_hmm_I), HMMOut_SNPs == AFFY_SNPs, ['hmm genotypes vs. Affy genotypes chrom ' num2str(chrom) ' Sample ' sample_names{cur_sample}]);
            end
        end

    end

end

% Calculate the error percents and the means
if(comp_hap_local)
    ErrorStruct.hap_local_common_snps = sum(ErrorStruct.hap_local_chr_common_snps_vec);
    ErrorStruct.hap_local_chr_precent_snps_vec = ErrorStruct.hap_local_chr_error_snps_vec ./ ...
        repmat(ErrorStruct.hap_local_chr_common_snps_vec, nsamples, 1); % Here starts a problem with dimensions! Need to fix
    ErrorStruct.hap_local_sample_precent_mean = sum(ErrorStruct.hap_local_chr_error_snps_vec,2)' ./ ErrorStruct.hap_local_common_snps;
    ErrorStruct.hap_local_chr_precent_mean = sum(ErrorStruct.hap_local_chr_error_snps_vec,1) ./ ErrorStruct.hap_local_chr_common_snps_vec;
    ErrorStruct.hap_local_percent_mean = sum(sum(ErrorStruct.hap_local_chr_error_snps_vec)) / (nsamples * ErrorStruct.hap_local_common_snps);
    if( length(ErrorStruct.hap_local_snp_ids)   < ErrorStruct.hap_local_common_snps)
        ErrorStruct.hap_local_snp_ids{ErrorStruct.hap_local_common_snps} = []; % Set the correct length
    end
end

if(comp_hap_hmm)
    ErrorStruct.hap_hmm_common_snps = sum(ErrorStruct.hap_hmm_chr_common_snps_vec);
    ErrorStruct.hap_hmm_chr_precent_snps_vec = ErrorStruct.hap_hmm_chr_error_snps_vec ./ ...
        repmat(ErrorStruct.hap_hmm_chr_common_snps_vec, nsamples, 1); % Here starts a problem with dimensions! Need to fix
    ErrorStruct.hap_hmm_sample_precent_mean = sum(ErrorStruct.hap_hmm_chr_error_snps_vec,2)' ./ ErrorStruct.hap_hmm_common_snps;
    ErrorStruct.hap_hmm_chr_precent_mean = sum(ErrorStruct.hap_hmm_chr_error_snps_vec,1) ./ ErrorStruct.hap_hmm_chr_common_snps_vec;
    ErrorStruct.hap_hmm_percent_mean = sum(sum(ErrorStruct.hap_hmm_chr_error_snps_vec)) / (nsamples * ErrorStruct.hap_hmm_common_snps);
end

if(comp_hap_affy)
    ErrorStruct.hap_affy_common_snps = sum(ErrorStruct.hap_affy_chr_common_snps_vec);
    ErrorStruct.hap_affy_chr_precent_snps_vec = ErrorStruct.hap_affy_chr_error_snps_vec ./ ...
        repmat(ErrorStruct.hap_affy_chr_common_snps_vec, nsamples, 1);
    ErrorStruct.hap_affy_sample_precent_mean = sum(ErrorStruct.hap_affy_chr_error_snps_vec,2)' ./ ErrorStruct.hap_affy_common_snps;
    ErrorStruct.hap_affy_chr_precent_mean = sum(ErrorStruct.hap_affy_chr_error_snps_vec,1) ./ ErrorStruct.hap_affy_chr_common_snps_vec;
    ErrorStruct.hap_affy_percent_mean = sum(sum(ErrorStruct.hap_affy_chr_error_snps_vec)) / (nsamples * ErrorStruct.hap_affy_common_snps);
end

if(comp_affy_hmm)
    ErrorStruct.affy_hmm_common_snps = sum(ErrorStruct.affy_hmm_chr_common_snps_vec);
    ErrorStruct.affy_hmm_chr_precent_snps_vec = ErrorStruct.affy_hmm_chr_error_snps_vec ./ ...
        repmat(ErrorStruct.affy_hmm_chr_common_snps_vec, nsamples, 1);
    ErrorStruct.affy_hmm_sample_precent_mean = sum(ErrorStruct.affy_hmm_chr_error_snps_vec,2)' ./ ErrorStruct.affy_hmm_common_snps;
    ErrorStruct.affy_hmm_chr_precent_mean = sum(ErrorStruct.affy_hmm_chr_error_snps_vec,1) ./ ErrorStruct.affy_hmm_chr_common_snps_vec;
    ErrorStruct.affy_hmm_percent_mean = sum(sum(ErrorStruct.affy_hmm_chr_error_snps_vec)) / (nsamples * ErrorStruct.affy_hmm_common_snps);
end