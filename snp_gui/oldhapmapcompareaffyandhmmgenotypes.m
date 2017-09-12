% This function reads the genotypes derived from the HAPMAP data, and
% compares it to the output of the HMM program, which is derived from the
% SNP-chip data that we analyzed. Clearly, one cannot expect to get 100%
% accuracy. A more reasonable goal is to get comparable, or better results,
% then those obtained by Affymetrix genotyping. To measure this, we also
% compare here the affymetrix genotypes with the HAPMAP genotypes, and
% obtain their accuracy
function ErrorStruct = HapMapCompareAffyAndHMMGenotypes(SampleNames, chip_type, user_dir, ...
    comp_hap_affy, comp_hap_hmm, comp_affy_hmm, use_strand_vec)
chroms = 1:22; % Currently do all the chromosomes
do_couples=1; joint_flag = 1; % Flags for getting marginals for HMM
nsamples = length(SampleNames);

if(comp_hap_hmm)
    ErrorStruct.hap_hmm_chr_common_snps_vec = zeros(1,22); % How many SNPs are on the chip and HAPMAP in each chrom
    ErrorStruct.hap_hmm_chr_error_snps_vec = zeros(nsamples,22); % For how many SNPs did we get a different genotype
end
if(comp_hap_affy)
    ErrorStruct.hap_affy_chr_common_snps_vec = zeros(1,22); % How many SNPs are on the chip and HAPMAP in each chrom for Affymetrix
    ErrorStruct.hap_affy_chr_error_snps_vec = zeros(nsamples,22); % For how many SNPs did we get a different genotype for Affymetrix
    ErrorStruct.affy_no_calls = zeros(nsamples,22);
end
if(comp_affy_hmm)
    ErrorStruct.affy_hmm_chr_common_snps_vec = zeros(1,22); % How many SNPs are on the chip and HMM in each chrom for Affymetrix
    ErrorStruct.affy_hmm_chr_error_snps_vec = zeros(nsamples,22); % For how many SNPs did we get a different genotype for Affymetrix
    ErrorStruct.affy_no_calls = zeros(nsamples,22);
end

% load the HAPMAP SNPs which appear on the chip (e.g. hind, xba or both) ...
if(comp_hap_affy || comp_hap_hmm)
    R=load('E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\CEU\HindXba_genotypes_chr_CEU_r21_nr_fwd.mat');
end

SampleInds = GetHapMapSampleIndex(SampleNames);
SampleBits = mod(SampleInds-1,30)+1;
SampleWords = ceil(SampleInds./30);

% Load Chip Annotations - We need the strand from here !!!
HindSNPChipAnnotStruct = load('..\..\Database\Hind_annot_data_hg17.mat', 'snp_ids', 'strand');

% Here we need to reverse the relevant SNPs according to strand
HindStrandSigns = zeros(1,length(HindSNPChipAnnotStruct.strand));
HindStrandSigns(strmatch('-', HindSNPChipAnnotStruct.strand)) = 1;

% Loop over all the samples
for cur_sample = 1:nsamples
    do_sample = cur_sample

    % Load the results of the Affymetrix genotypes
    if(comp_hap_affy || comp_affy_hmm)
        %%    eval(['Aff = load(''' user_dir SampleNames{cur_sample} '_hind.mat'', ''snp_id_hind'', ''genotype_vec_hind'');']);
        eval(['Aff = load(''' user_dir SampleNames{cur_sample} '_hind.mat'');']);  % Should we load also the allele ratios? Yes if we want to plot stuff!
    end

    % Load the results of the HMM
    if(comp_hap_hmm || comp_affy_hmm)
        eval(['load(''' user_dir 'hmm_out\' SampleNames{cur_sample} '_hind_hmm_out.mat'');']);  %% '_hind_UNIPLACE_hmm_out.mat'
    end


    % Set the relevant string and strand - NOTE: We do it for the HAPMAP data and NOT for the affy and/or HMM genotypes !!!
    [snp_ids_intersect I J] = intersect(HindSNPChipAnnotStruct.snp_ids, Aff.snp_id_hind);
    %%     eval(['Aff.genotype_vec_' chip_type '(J(find(HindStrandSigns(I)))) = 4-Aff.genotype_vec_' chip_type '(J(find(HindStrandSigns(I))));']);
    %%     eval(['Aff.genotype_vec_' chip_type '(find(Aff.genotype_vec_' chip_type ' == 0)) = 4;']);

    %% The same should be done for the Xba later ...
    %%    XbaSNPChipAnnotStruct = load('..\..\Database\Xba_annot_data_hg17.mat', 'snp_ids', 'strand');
    %% XbaHT_SIGNS = zeros(1,length(XbaSNPChipAnnotStruct.snp_ids));
    %% XbaHT_SIGNS(strmatch('-', XbaSNPChipAnnotStruct.strand)) = 1;

    % Go chrom-by-chrom as to not take too much memory
    ErrorStruct.hap_affy_ErrorMat{cur_sample} = zeros(4);
    ErrorStruct.hap_hmm_ErrorMat{cur_sample} = zeros(4);
    ErrorStruct.affy_hmm_ErrorMat{cur_sample} = zeros(4);

    for chrom=chroms
        do_chrom = chrom

        if(comp_hap_hmm || comp_affy_hmm) % isolate to save time ...
            [VitStruct ProbsStruct] = ...
                GetBestMarginalPredictions(HMMOutStruct.SNPsProbsMat{chrom}, HMMOutStruct.ModelParams, do_couples, joint_flag);
            %%            [VitStruct2 ProbsStruct2] = ...
            %%                GetBestMarginalPredictions(HMMOutStruct.SNPsProbsMat{chrom}, HMMOutStruct.ModelParams, 0, joint_flag); % Try without couples ..
        end
        if(comp_hap_hmm)
            [snp_ids_intersect I J] = intersect(HMMOutStruct.SNPsIDs{chrom}, R.SnpsIDs{chrom}); % Get the common SNPs
            ErrorStruct.hap_hmm_chr_common_snps_vec(chrom) = length(I);


            HAPMAP_SNPs = 2*bitget(R.SnpsData{chrom}(J,SampleWords(cur_sample)), SampleBits(cur_sample)) + ...
                bitget(R.SnpsData{chrom}(J,SampleWords(cur_sample)+3), SampleBits(cur_sample)); % Take the first patient from HAPMAP. Bad Calls are considered as '-1'
            HMMOut_SNPs = VitStruct.joint_genotype(I);  % Take the patient from the HMM Output ..
            %% Correct for the strand ... should we ?? NO %[strand_snp_ids_intersect strand_I strand_J] = intersect(HindSNPChipAnnotStruct.snp_ids, HMMOutStruct.SNPsIDs{chrom}); % Get the common SNPs
            %%  [snp_ids_intersect strand_I strand_J] = intersect(HindSNPChipAnnotStruct.snp_ids, HMMOutStruct.SNPsIDs{chrom}(I));
            %%  HMMOut_SNPs(strand_J(find(HindStrandSigns(strand_I)))) = 3-HMMOut_SNPs(strand_J(find(HindStrandSigns(strand_I))));
            HAPMAP_SNPs(find(bitget(R.SnpsBadCalls{chrom}(J,1), 1))) = -1; % Bad calls NN SNPs are set to minus one

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
            [hap_affy_snp_ids_intersect hap_affy_I hap_affy_J] = intersect(Aff.snp_id_hind, R.SnpsIDs{chrom}); % Get the common SNPs
            ErrorStruct.hap_affy_chr_common_snps_vec(chrom) = length(hap_affy_I);

            HAPMAP_SNPs = 2*bitget(R.SnpsData{chrom}(hap_affy_J,SampleWords(cur_sample)), SampleBits(cur_sample)) + ...
                bitget(R.SnpsData{chrom}(hap_affy_J,SampleWords(cur_sample)+3), SampleBits(cur_sample)); % Take the first patient from HAPMAP
            if(use_strand_vec(2))             % Correct for reverse strand
                [strand_snp_ids_intersect strand_I strand_J] = intersect(HindSNPChipAnnotStruct.snp_ids,  R.SnpsIDs{chrom}(hap_affy_J)); % Get the common SNPs
                HAPMAP_SNPs(strand_J(find(HindStrandSigns(strand_I)))) = 3-HAPMAP_SNPs(strand_J(find(HindStrandSigns(strand_I))));
            end
            HAPMAP_SNPs(HAPMAP_SNPs ==2) = 1;
            HAPMAP_SNPs(find(bitget(R.SnpsBadCalls{chrom}(hap_affy_J,1), 1))) = -1; % Bad calls NN SNPs are set to minus one

            AFFY_SNPs = Aff.genotype_vec_hind(hap_affy_I)-1; % Get the SNPs that were obtained by affymetrix software from their chips ...
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
            [affy_hmm_snp_ids_intersect affy_hmm_I affy_hmm_J] = intersect(Aff.snp_id_hind, HMMOutStruct.SNPsIDs{chrom}); % Get the common SNPs
            ErrorStruct.affy_hmm_chr_common_snps_vec(chrom) = length(affy_hmm_I);
            AFFY_SNPs = Aff.genotype_vec_hind(affy_hmm_I)-1; % Get the SNPs that were obtained by affymetrix software from their chips ...
            AFFY_SNPs(AFFY_SNPs == 3) = -1; % The No-Calls
            AFFY_SNPs(AFFY_SNPs == 2) = 3; % The BB's
            HMMOut_SNPs = VitStruct.joint_genotype(affy_hmm_J);  % Take the patient from the HMM Output ..
            HMMOut_SNPs(HMMOut_SNPs == 2) = 1; % AB and BA are the same ..

            % Correct for reverse strand - should we ? YES
            if(use_strand_vec(3))
                [snp_ids_intersect strand_I strand_J] = intersect(HindSNPChipAnnotStruct.snp_ids, HMMOutStruct.SNPsIDs{chrom}(affy_hmm_J));
                HMMOut_SNPs(strand_J(find(HindStrandSigns(strand_I)))) = 3-HMMOut_SNPs(strand_J(find(HindStrandSigns(strand_I))));
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
                PlotAlleleRatios(Aff.copy_num_vec_hind(hap_affy_I), Aff.allele_ratio_vec_hind(hap_affy_I), max(0,HAPMAP_SNPs), ['HAPMAP genotypes chrom ' num2str(chrom) ' Sample ' SampleNames{cur_sample}]);
                PlotAlleleRatios(Aff.copy_num_vec_hind(hap_affy_I), Aff.allele_ratio_vec_hind(hap_affy_I), Aff.genotype_vec_hind(hap_affy_I), ['affy genotypes chrom ' num2str(chrom) ' Sample ' SampleNames{cur_sample}]);
                PlotAlleleRatios(Aff.copy_num_vec_hind(hap_affy_I), Aff.allele_ratio_vec_hind(hap_affy_I), AFFY_SNPs == HAPMAP_SNPs, ['affy vs. hap genotypes  Sample ' SampleNames{cur_sample}]);
            end

            if(comp_affy_hmm || comp_hap_affy)
                PlotAlleleRatios(Aff.copy_num_vec_hind, Aff.allele_ratio_vec_hind, Aff.genotype_vec_hind, ['affy genotypes  Sample ' SampleNames{cur_sample}]);
            end
           if(comp_hap_hmm)
               PlotAlleleRatios(Aff.copy_num_vec_hind(affy_hmm_I), Aff.allele_ratio_vec_hind(affy_hmm_I), HMMOut_SNPs == HAPMAP_SNPs, ['hmm genotypes vs. HapMap genotypes chrom ' num2str(chrom) ' Sample ' SampleNames{cur_sample}]);
           end
            
            if(comp_affy_hmm)
                PlotAlleleRatios(Aff.copy_num_vec_hind(affy_hmm_I), Aff.allele_ratio_vec_hind(affy_hmm_I), max(0, AFFY_SNPs), ['affy genotypes chrom ' num2str(chrom) ' Sample ' SampleNames{cur_sample}]); % Aff.genotype_vec_hind(affy_hmm_I)
                PlotAlleleRatios(Aff.copy_num_vec_hind(affy_hmm_I), Aff.allele_ratio_vec_hind(affy_hmm_I), HMMOut_SNPs, ['hmm genotypes chrom ' num2str(chrom) ' Sample ' SampleNames{cur_sample}]);
                PlotAlleleRatios(Aff.copy_num_vec_hind(affy_hmm_I), Aff.allele_ratio_vec_hind(affy_hmm_I), HMMOut_SNPs == AFFY_SNPs, ['hmm genotypes vs. Affy genotypes chrom ' num2str(chrom) ' Sample ' SampleNames{cur_sample}]);
            end
        end

    end

end

% Calculate the error precents and the means
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
