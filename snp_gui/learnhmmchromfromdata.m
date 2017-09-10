% Written by Or Zuk 6/2007
%
% A function for learning the HMM models from the data. This function does most
% of the work when calling 'Run HMM' from the SNP GUI. The function saves
% all the output in files - one file per sample for display and one hmm_out
% file.
%
% NEW: We try the SNP-specific HMM here. No learning at all
%
% The inputs:
% user_dir - working directory with input files
% SampleNames - names of samples
% LDStruct - structure containing Linkeage-Disequilibrium statistics from Hapmap
% SNPChipAnnotStruct - Annotations for the specific SNP chip
% HMMParamsStruct - parameters for the HMM running
%
% The outputs:
% OutputFilesNames - the names of the files in which the function saves the HMM output.
function OutputFilesNames = ...
    LearnHMMChromFromData(user_dir, SampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct)

ttt_all = cputime;
ttt = cputime;
fid = fopen('time_debug.txt','wt');

% Various data types
AssignAllGlobalConstants();
HMMParamsStruct.use_affy_genotypes = 0; % Should be eliminated - input as a button ...
HMMParamsStruct.SmoothWindow = 30; % Size of window to do smoothing with
hmm_out_dir = fullfile(user_dir, '\hmm_out');
find_path_Viterbi=1; use_locations=0; do_fold_change=0; do_determine_x_dim=0; derich = 1.0/120.0; % Relative Derichlet correction
chip_type = lower(SNPChipAnnotStruct.chip);
num_samples = length(SampleNames);
SPECIAL_MODELS_FLAG = 1; % Must have the special SNPs model
RelevantInds = [1:6 17:22 33:38 49:54 65:70 81:86]; % This is good for the case of 3 levels (amp/normal/deleted)
do_couples = 0; % 1 take into acount adjacent couples in analysis
joint_flag=1; % joint flag - what is this???
plot_for_debugging=0; % enable/disable several plots
HMMParamsStruct.SNP_specific = 1; % NEW! We have different parameters for different SNPs
% load the correction for the data at once
SHIFTS = load(fullfile(user_dir, ['AllSamplesMat_' chip_type '.mat']), 'NormalMeanIntensities', 'SampleNames');
SHIFTS.NormalMultIntensities = 2 ./ SHIFTS.NormalMeanIntensities; % this is the multiplicative correction
[SHIFTS.SampleNames I J] = intersect_order_by_first_gr(SampleNames, SHIFTS.SampleNames);
SHIFTS.NormalMultIntensities = SHIFTS.NormalMultIntensities(J);

r_mat = zeros(4);
%% r_mat(2,:) = [ 0.6886    1.0000    1.1103    1.2446];
r_mat(2,:) = [ 0.6886    1.0000    1.2446    1.5000];

% Load the RLMM parameters
if(HMMParamsStruct.SNP_specific)
    HMMParamsStruct.learn_model_EM = 0; % Here we don't learn the model but get it ...
    load(fullfile('../database', ['RLMM_' pop_str_vec{HMMParamsStruct.hapmap_population} '_' chip_type '_' genome_assembly '.mat']));
    RLMM = RLMM_ComputeAuxillaryParams(RLMM); % get all parameters - just to be sure ...
end

ttt = cputime - ttt
fprintf(fid,'start time %lf\n',ttt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Old - first loop on samples - probably not neccessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ALL_MAT = {};
% % % for(cur_sample = 1:num_samples) % Outer loop on samples ....
% % %     sample_name = SampleNames{cur_sample}
% % %     ttt = cputime;
% % %     HMMOutStruct = {};
% % %     for i=1:num_chr_to_run
% % %         HMMOutStruct.SNPsIDs{i} = {};
% % %     end
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     % Just make all the indexology right here
% % %     [TMP_ALL_MAT CHROM_MATS HMMOutStruct] = LoadAndIntersectData(user_dir, sample_name, LDStruct, SNPChipAnnotStruct, ...
% % %         HMMParamsStruct, HMMOutStruct);
% % %     % New! Correct data to be around 2
% % %     TMP_ALL_MAT.data = TMP_ALL_MAT.data .*  SHIFTS.NormalMultIntensities(cur_sample);
% % %     TMP_ALL_MAT.data_A = TMP_ALL_MAT.data_A .*  SHIFTS.NormalMultIntensities(cur_sample);
% % %     TMP_ALL_MAT.data_B = TMP_ALL_MAT.data_B .*  SHIFTS.NormalMultIntensities(cur_sample);
% % %     for i=1:length(CHROM_MATS.data_A)
% % %         CHROM_MATS.data_A{i} = CHROM_MATS.data_A{i} .* SHIFTS.NormalMultIntensities(cur_sample);
% % %         CHROM_MATS.data_B{i} = CHROM_MATS.data_B{i} .* SHIFTS.NormalMultIntensities(cur_sample);
% % %         HMMOutStruct.SNPsIDs{i} = {};
% % %     end
% % % 
% % %     num_snps = length(TMP_ALL_MAT.loc);
% % % 
% % %     if(cur_sample==1) % On the first sample prepare the matrices
% % %         ALL_MAT.data_genotype_A = zeros(num_samples, num_snps);
% % %         ALL_MAT.data_genotype_B = zeros(num_samples, num_snps);
% % %         ALL_MAT.data_genotype_AB = zeros(num_samples, num_snps);
% % % 
% % %         ALL_MAT.snp_ids = TMP_ALL_MAT.snp_ids;
% % %         ALL_MAT.loc = TMP_ALL_MAT.loc;
% % %         ALL_MAT.strand = TMP_ALL_MAT.strand;
% % %         ALL_MAT.data_A = zeros(num_samples, num_snps);
% % %         ALL_MAT.data_B = zeros(num_samples, num_snps);
% % %     end
% % %     % Copy the intensities
% % %     ALL_MAT.data_A(cur_sample,:) = TMP_ALL_MAT.data_A;
% % %     ALL_MAT.data_B(cur_sample,:) = TMP_ALL_MAT.data_B;
% % % 
% % % 
% % % 
% % %     HMM_MODEL = {};
% % % 
% % %     ttt = cputime - ttt
% % %     fprintf(fid,'first loop sample %ld  time %lf\n',cur_sample, ttt);
% % % 
% % % end % loop for fast PHASING (and currently also classifying ..)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New: Classify using RLMM:                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(fullfile('../database', ['RLMM_' pop_str_vec{HMMParamsStruct.hapmap_population} '_' chip_type '_' genome_assembly '.mat']));
% [ALL_MAT.data_genotype_AB Mahalahonis] = ...
%     RLMM_ClassifyUsingGaussianParams(ALL_MAT.data_A', ALL_MAT.data_B', RLMM.MuMats, RLMM.SigmaMats);
% PlotAlleleRatios(ALL_MAT.data_A(1,:) + ALL_MAT.data_B(1,:), ALL_MAT.data_A(1,:) ./ ALL_MAT.data_B(1,:), ...
%      ALL_MAT.data_genotype_AB(:,1) , 'One Sample RLMM Out');
% PlotAlleleRatios(ALL_MAT.data_A(1,1:1000) + ALL_MAT.data_B(1,1:1000), ALL_MAT.data_A(1,1:1000) ./ ALL_MAT.data_B(1,1:1000), ...
%      ALL_MAT.data_genotype_AB(1:1000,1) , 'One Sample RLMM Out 1000 SNPs');
% ErrorStructRLMM = HapMapCompareAffyAndHMMGenotypes(SampleNames, chip_type, [user_dir '\'], ALL_MAT, 0, 1, 0, [0 1 0 0], ...
%     HMMParamsStruct.hapmap_dir, HMMParamsStruct.hapmap_population, HMMParamsStruct.hapmap_version);
% save(['BLRMM_Hapmap' chip_type pop_str_vec{HMMParamsStruct.hapmap_population} 'ErrorStructXba.mat'], 'ErrorStructRLMM');
% TTT = 999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Perform PCA special for the genotype
plot_pca=0;
if((num_samples > 1) && (plot_pca == 1) && HMMParamsStruct.use_affy_genotypes)% No point in doing PCA for one sample. -1 means all chromosomes
    DataLabels = SampleNames;
    for j=1:length(SampleNames)
        DataLabels{j} = SampleNames{j}(end);
    end
    DisplaySamplesPCA(ALL_MAT, 1, DataLabels, -1, 1000);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform PHASING by using an outside program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FastPhase_dir = fullfile(user_dir, '\FastPhase');
do_phase=0;
if(do_phase)
    if(~exist(FastPhase_dir)) % Make the directory if needed
        eval(['mkdir ' FastPhase_dir]);
    end
    FastPhaseInputFile = fullfile(user_dir, 'FastPhase', 'AllSamples_FastPhaseInput.txt');
    FastPhaseOutputFile = fullfile(user_dir, 'FastPhase', 'AllSamples');
    PrepareInputFileForFastPhase(ALL_MAT, FastPhaseInputFile, SampleNames); % Prepare input for Stephens' fast Phasing algorithm
    cur_dir = pwd; eval(['cd ' user_dir '\FastPhase']); % Go to user directory
    FastPhaseOutputFile = [FastPhaseOutputFile '_hapguess_switch.out'];
    if(~exist(FastPhaseOutputFile)) % Run Stephens' fast Phasing algorithm if needed - this is not fast at all!!!
        system([cur_dir '\fastPHASE -T1 -C10 -oAllSamples AllSamples_FastPhaseInput.txt']);
    end
    eval(['cd ' cur_dir]); % Go back to our directory
end
%% M=ReadOutputFileFromFastPhase(FastPhaseOutputFile, SampleNames); % Read output of Stephens' fast Phasing algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of PHASING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

average_copy_num_mat = zeros(0, num_samples, 'single');
copy_num_mat = zeros(0, num_samples, 'single');
genotype_mat = zeros(0, num_samples, 'single');
data_snp_ids = cell(0,1);

% load the matrix withh all the data for all the samples
OUTLIERS = load(fullfile(user_dir, ['AllSamplesMat_' chip_type '.mat']), 'array_outlier_vec', 'SampleNames');
[OUT_INT OUT_I OUT_J] = intersect_order_by_first_gr(SampleNames, OUTLIERS.SampleNames); 
array_outlier_vec = OUTLIERS.array_outlier_vec(OUT_J);
beta = 10000; p_vec = double(max(0.01 .* exp(-beta.*array_outlier_vec), MIN_P_TRANSITION_VAL));
%% p_vec = 0.01 .* ones(1,num_samples); 

for cur_sample = 1:num_samples % Outer loop on samples ....
    sample_name = SampleNames{cur_sample};
    ttt = cputime;
    HMMOutStruct = {};
    hmm_out_file_name = fullfile('hmm_out',[sample_name '_' chip_type '_hmm_out']);
    for i=1:24
        HMMOutStruct.SNPsIDs{i} = {};
    end
    if( mod(cur_sample,10) == 0 )
        TTTT = 919191;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ALL_MAT CHROM_MATS HMMOutStruct] = LoadAndIntersectData(user_dir, sample_name, LDStruct, SNPChipAnnotStruct, ...
        HMMParamsStruct, HMMOutStruct);
    % New! Correct data to be around 2
    ALL_MAT.data = ALL_MAT.data .*  SHIFTS.NormalMultIntensities(cur_sample);
    ALL_MAT.data_A = ALL_MAT.data_A .*  SHIFTS.NormalMultIntensities(cur_sample);
    ALL_MAT.data_B = ALL_MAT.data_B .*  SHIFTS.NormalMultIntensities(cur_sample);
    for i=1:length(CHROM_MATS.data_A)
        CHROM_MATS.data_A{i} = CHROM_MATS.data_A{i} .* SHIFTS.NormalMultIntensities(cur_sample);
        CHROM_MATS.data_B{i} = CHROM_MATS.data_B{i} .* SHIFTS.NormalMultIntensities(cur_sample);
        HMMOutStruct.SNPsIDs{i} = {};
    end

    HMM_MODEL = {}; %% HMM_MODEL.PLACE_M3 = PLACE_M3;
    %%    HMM_MODEL.PLACE_M3(find(ALL_MAT.strand==1),:) = 1-HMM_MODEL.PLACE_M3(find(ALL_MAT.strand==1),:);

    % Always update the PLACE_M according to the whole data. Cond. prob.: PLACE_M(i,j) = Pr(S_{n+1} = j | S_n = i)
    HMMParamsStruct.derich = derich;
    HMM_MODEL = LDToPlaceM(LDStruct, HMM_MODEL, HMMParamsStruct); % update the PLACE_M matrix according to Linkage-Disequilibrium


    %%     % Compute mutual information for place matrices
    %%     for i=1:num_chr_to_run
    %%         LDStruct.LD{i}.MutualInfo = -entropy(LDStruct.LD{i}.PairMat') + ...
    %%             entropy([LDStruct.LD{i}.PairMat(:,1) + LDStruct.LD{i}.PairMat(:,2), ...
    %%             LDStruct.LD{i}.PairMat(:,3) + LDStruct.LD{i}.PairMat(:,4)]') + ...
    %%             entropy([LDStruct.LD{i}.PairMat(:,1) + LDStruct.LD{i}.PairMat(:,3), ...
    %%             LDStruct.LD{i}.PairMat(:,2) + LDStruct.LD{i}.PairMat(:,4)]');
    %%         LDStruct.LD{i}.MutualInfo = max(LDStruct.LD{i}.MutualInfo, 0); % Avoid very small negative numbers
    %%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get some neccessary model parameters
    HMM_MODEL.PLACE_FLAG = 1;
    HMM_MODEL.SPECIAL_MODELS_FLAG = 1;
    % Note : We assume a stationary distribution here so PI is determined accordingly
    % Learn Again ! (This should be unneccessary) - we learn a different
    % model for each sample, but the model is the SAME for all the chromosomes !!!!
    HMM_MODEL.PLACE_M_C = HMM_MODEL.PLACE_M;
    HMM_MODEL.PLACE_M_C(:,1) = HMM_MODEL.PLACE_M(:,2);
    HMM_MODEL.PLACE_M_C(:,2) = HMM_MODEL.PLACE_M(:,1);
    HMM_MODEL.PLACE_M_C(:,3) = HMM_MODEL.PLACE_M(:,4);
    HMM_MODEL.PLACE_M_C(:,4) = HMM_MODEL.PLACE_M(:,3);
    % Switch 2 and 3
    HMM_MODEL.PLACE_M_C(:,2) = HMM_MODEL.PLACE_M(:,4);
    HMM_MODEL.PLACE_M_C(:,3) = HMM_MODEL.PLACE_M(:,1);
    ttt_EM = cputime;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Learn a model for each chromosomes of each patient
    % here learn the model parameters
    if(HMMParamsStruct.learn_model_EM)
        % Dummy variables for calling c-function
        mean_vec_rep = zeros(num_samples,HMMParamsStruct.x_dim);
        std_vec_rep = zeros(num_samples,HMMParamsStruct.x_dim);
        % Note : Here we try to learn the SIZE of the model (the number of states X can attain !!) Currently not used!
        if(do_determine_x_dim)
            max_x_dim = HMMParamsStruct.x_dim; % Now this is the maximal value, while the actual value might be smaller
            [K_min, pen_dist_vec, params_cell, logscores] = ...
                HMM_K_by_penalized_dist_EM(HMM_chrom_data, 1:max_x_dim, HMMParamsStruct.num_EM_iters, HMMParamsStruct.num_EM_starting_points, HMMParamsStruct.EM_tolerance);

            % Now copy all the models !!!
            for Ksize = 1:max_x_dim
                HMM_MODEL.M =  params_cell{Ksize,1};
                HMM_MODEL.N =  ones(1, Ksize); % Currently we have no mixtures !!!!
                HMM_MODEL.MEW =  params_cell{Ksize,2}(:,1);
                HMM_MODEL.SIGMA =  params_cell{Ksize,2}(:,2);
                HMM_MODEL.PI = params_cell{Ksize,3};
                HMM_MODEL.LogScore = logscores(Ksize);  % Now save also the output score !!!
            end % End loop on Ksize
            % Tried moving it a bit upwards
            HMM_MODEL.X_DIM = Ksize;

            %  HMMParamsStruct.x_dim = K_min;
        else  % Here we know the model size
            K_min = HMMParamsStruct.x_dim; % HMMParamsStruct.x_dim
            Ksize = HMMParamsStruct.x_dim;
        end  % if determine_x_dim

        [HMM_MODEL.PI HMM_MODEL.M HMM_MODEL.N HMM_MODEL.MEW ...
            HMM_MODEL.SIGMA HMM_MODEL.LogScore] = ...
            TrainHMMFromDataEMMatlab(ALL_MAT.data_A, ALL_MAT.data_B, ALL_MAT.loc, ...
            use_locations, K_min, HMMParamsStruct.y_dim, ...
            1, do_fold_change, mean_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M', ...
            SPECIAL_MODELS_FLAG, eye(Ksize), 0, ...
            HMMParamsStruct.num_EM_iters, HMMParamsStruct.num_EM_starting_points, HMMParamsStruct.EM_tolerance); % Here we put the upper-bounds. Currently they're not used !!!

        % set parameters instead of training
        %         HMM_MODEL.PI = ones(3,1)./3;
        %         HMM_MODEL.M = [0.9504 0.0283 0.0213; 0.0131 0.9565 0.0304; 0.0376 0.0624 0.9000];
        %         HMM_MODEL.N = ones(5,1);
        %         HMM_MODEL.MEW = [0.1 0.8 2.0 2.5 3.0]';
        %         HMM_MODEL.SIGMA = [0.5 0.8 1.0 1.5 2.0]';
        %         HMM_MODEL.LogScore = -100;

        ttt_EM = cputime - ttt_EM
        CurrentLogScore = HMM_MODEL.LogScore
    else % Do not learn parameters
        K_min = HMMParamsStruct.x_dim; % HMMParamsStruct.x_dim
        Ksize = HMMParamsStruct.x_dim;
        HMM_MODEL.PI = ones(HMMParamsStruct.x_dim,1)./HMMParamsStruct.x_dim; pp = 0.9;
        HMM_MODEL.M = pp * eye(HMMParamsStruct.x_dim) + ((1-pp)/(HMMParamsStruct.x_dim-1)) * (1- eye(HMMParamsStruct.x_dim));
        HMM_MODEL.N = ones(2*HMMParamsStruct.x_dim-1,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % These are not used anyway ...                      %
        HMM_MODEL.MEW = [0.1 0.8 2.0 2.5 3.0]';              %
        HMM_MODEL.SIGMA = [0.5 0.8 1.0 1.5 2.0]';            %
        HMM_MODEL.LogScore = -99;                            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Need to check that we're not missing anything here !!!
        %%   load(fullfile(user_dir, hmm_out_file_name)); % we assume that this exists ...
        %% HMM_MODEL = HMMOutStruct.ModelParams; % Take the already saved model (why?)
        HMM_MODEL = LDToPlaceM(LDStruct, HMM_MODEL, HMMParamsStruct); % update the PLACE_M matrix according to Linkage-Disequilibrium

    end   % if HMMParamsStruct.learn_model_EM - finished learning part

    ttt_VIT = cputime;
    chrom_start_ind=1;

    ttt = cputime - ttt
    fprintf(fid,'sample %ld before Viterby time %lf\n',cur_sample, ttt);


    p = p_vec(cur_sample);  HMM_MODEL.M2 = p .* HMM_MODEL.M + (1-p) .* eye(HMMParamsStruct.x_dim); % take a convex combination, try a close to 'non-switching' matrix

    % Here do the Viterbi and Forward/Backward part
    for i=HMMParamsStruct.ChromosomesToRun %%% num_chroms % go over all chromosomes
        NowDoViterbiChromosome = i
        ttt = cputime;
        if(~isempty(LDStruct.LD{i}.RS_IDs))
            num_genes_in_chrom(i) = length(LDStruct.LD{i}.RS_IDs);

            HMM_MODEL_M2_IS = HMM_MODEL.M2     
            
            % Problem: We may have some (very few) SNPs which are not in the HAPMAP but on the chip.
            % We want to insert them also so we need to 'strech' the M's (or 'remove' the data).
            % This actualy does the 'classification'
            if(find_path_Viterbi)
                Viterbi_Path{i} = zeros(num_samples, num_genes_in_chrom(i) );  % how many genes are in the current chromosome
                mean_vec = zeros(1,num_genes_in_chrom(i)); std_vec = 1+zeros(1,num_genes_in_chrom(i));  % Dummy vectors again ...
                HMM_MODEL.Y_TYPE = GAUSSIAN;
                HMM_MODEL.y_dim = 1; % No mixture, only one gaussian ...
                HMM_MODEL.x_dim = HMMParamsStruct.x_dim;
%%                HMM_MODEL.M2 = p .* HMM_MODEL.M + (1-p) .* eye(HMMParamsStruct.x_dim); % take a convex combination, try a close to 'non-switching' matrix
%%                p = 0.01; HMM_MODEL.M2 = p .* HMM_MODEL.M + (1-p) .* eye(HMMParamsStruct.x_dim); % take a convex combination, try a close to 'non-switching' matrix

                chrom_end_ind = chrom_start_ind + length(CHROM_MATS.data_A{i})-1;

                start_ind = 1; end_ind = length(CHROM_MATS.data_A{i});
                L = CHROM_MATS.loc{i}(start_ind:end_ind);
                use_x_flag = 1; % This flag says if we want to calculate log(Pr(X,Y)) or just log(Pr(Y))
                %             DataLogLikelihood2 = ...
                %                 ComputeHMMLogLikelihood(Viterbi_Path{i}, CHROM_MATS.data_A{i}, CHROM_MATS.data_B{i}, use_x_flag, L, ...
                %                 use_locations, HMM_MODEL.PI, HMM_MODEL.M2, HMM_MODEL.N, ...
                %                 HMM_MODEL.MEW', HMM_MODEL.SIGMA, HMM_MODEL.PLACE_FLAG, ...
                %                 HMM_MODEL.SPECIAL_MODELS_FLAG, mean_vec, std_vec, ...
                %                  SHORT_PLACE_M_C');

                SHORT_PLACE_M = HMM_MODEL.PLACE_M(chrom_start_ind:chrom_end_ind,:); chrom_start_ind=chrom_end_ind+1;
                SHORT_PLACE_M_C = SHORT_PLACE_M; % Adjust dimension
                SHORT_PLACE_M_C(:,1) = SHORT_PLACE_M(:,2);
                SHORT_PLACE_M_C(:,2) = SHORT_PLACE_M(:,1);
                SHORT_PLACE_M_C(:,3) = SHORT_PLACE_M(:,4);
                SHORT_PLACE_M_C(:,4) = SHORT_PLACE_M(:,3);
                % Switch 2 and 3
                SHORT_PLACE_M_C(:,2) = SHORT_PLACE_M(:,4);
                SHORT_PLACE_M_C(:,3) = SHORT_PLACE_M(:,1);

                % NEW: Here we may choose to use the SNP-specific functions ..
                HMM_MODEL.SNP_specific = HMMParamsStruct.SNP_specific;

                if(HMMParamsStruct.SNP_specific)
                    % perform intersection with RLMM to get the correct SNPs moments
                    [mean_mat std_mat std_inv_mat] = ExpandSNPsMoments(RLMM, CHROM_MATS.snp_ids{i}, r_mat);
                    % start with the Viterby - Note: All is Double !!!
                    [Viterbi_Path{i} HMMOutStruct.SNPsProbsMat{i}] = ...
                        FindBestPathViterbi(double(CHROM_MATS.data_A{i}), double(CHROM_MATS.data_B{i}), double(L), ...
                        use_locations, HMM_MODEL.PI, HMM_MODEL.M2, HMM_MODEL.N, ...
                        HMM_MODEL.MEW', HMM_MODEL.SIGMA, HMM_MODEL.PLACE_FLAG, ...
                        HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.SNP_specific, double(mean_mat), double(std_mat), double(std_inv_mat), ...
                        SHORT_PLACE_M_C');

                else
                    [Viterbi_Path{i} HMMOutStruct.SNPsProbsMat{i}] = ...
                        FindBestPathViterbi(CHROM_MATS.data_A{i}, CHROM_MATS.data_B{i}, L, ...
                        use_locations, HMM_MODEL.PI, HMM_MODEL.M2, HMM_MODEL.N, ...
                        HMM_MODEL.MEW', HMM_MODEL.SIGMA, HMM_MODEL.PLACE_FLAG, ...
                        HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.SNP_specific, mean_vec, std_vec, std_inv_vec, ...
                        SHORT_PLACE_M_C');
                end
                [VitStruct{i} ProbsStruct{i}] = ...
                    GetBestMarginalPredictions(HMMOutStruct.SNPsProbsMat{i}, HMM_MODEL, do_couples, joint_flag);
                clear Viterbi_Path{i}; % Save some more memory ..
            end  % Viterbi
        end % check if empty chromosome
        ttt = cputime - ttt
        fprintf(fid,'sample %ld  Viterby chrom %ld time %lf\n',cur_sample, i, ttt);


    end  % loop on chromosomes (i)
    ttt = cputime;
    ttt_VIT = cputime - ttt_VIT
    clear R;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Plot for debugging: (should be removed later)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    if(plot_for_debugging)
        ttt_PLOT = cputime;
        PlotForDebugging(CHROM_MATS, user_dir, sample_name, chip_type, ...
            SNPChipAnnotStruct, HMMParamsStruct, HMM_MODEL, VitStruct, ProbsStruct, ALL_MAT);
        ttt_PLOT = cputime - ttt_PLOT
    end

    ttt_MIX = cputime;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Perform simple MOG classification
    %     INIT_S_2D{1} = [1 0; 0 1]; INIT_S_2D{2} = [1 0; 0 1];  INIT_S_2D{3} = [1 0; 0 1];
    %     INIT_M_2D = [1 0; 1 1; 0 1]; INIT_P_2D = ones(1,3)/3;
    %     [P_2D,M_2D,SIGMA_2D, LogLike_2d]= ...
    %         MixtureOfGaussiansMultiDimGivenInit([TMP_ALL_MAT.data_A' TMP_ALL_MAT.data_B']', ...
    %         3,25, INIT_P_2D, INIT_M_2D, INIT_S_2D);
    %     C = MixtureOfGaussiansMultiDimClassify([TMP_ALL_MAT.data_A' TMP_ALL_MAT.data_B']', P_2D, M_2D, SIGMA_2D); C=C-1; C(C==2)=3;
    %     PlotAlleleRatios(TMP_ALL_MAT.data_A + TMP_ALL_MAT.data_B, ...
    %         TMP_ALL_MAT.data_A ./ TMP_ALL_MAT.data_B, C, 'MY MOG genotypes ');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    if(isfield(HMM_MODEL, 'PLACE_M_C')) % try to save some memory
        HMM_MODEL = rmfield(HMM_MODEL, 'PLACE_M_C');
    end
    for i=HMMParamsStruct.ChromosomesToRun %%% num_chroms % go over all chromosomes
        if(~isempty(LDStruct.LD{i}.RS_IDs))
            HMMOutStruct.SNPsProbsMat{i} = sparse(HMMOutStruct.SNPsProbsMat{i} .* (HMMOutStruct.SNPsProbsMat{i} > EPSILON));
        else
            HMMOutStruct.SNPsProbsMat{i} = [];
            VitStruct{i} = [];
            ProbsStruct{i} = [];
        end
    end
    [display_file_name DispStruct] = DisplayFormatSave(user_dir, sample_name, chip_type, genome_assembly, CHROM_MATS, ...
        HMMParamsStruct, VitStruct, ProbsStruct);
    ttt = cputime - ttt
    fprintf(fid,'sample %ld  display save time %lf\n',cur_sample, ttt);
    ttt = cputime;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_chr = length(DispStruct.Chrom);
    if(cur_sample==1)
        chr_num_snps_vec = zeros(1, num_chr);
        num_of_blank_SNPs = zeros(1, num_chr);
    end
    snp_ind=1;
    for j = 1:num_chr
        if(~isempty(DispStruct.Chrom{j}))
            if(cur_sample==1)
                % here should be changed by OR
                data_snp_ids = concat_cells(data_snp_ids, DispStruct.Chrom{j}.SNPsIDs', 2);
                chr_num_snps_vec(j) = size(DispStruct.Chrom{j}.SNPsIDs, 2);
                num_of_blank_SNPs(j)=length(find(strcmp(deblank(data_snp_ids),'')))-sum(num_of_blank_SNPs); %count number of blank SNPs 21/03/2007
            end
            genotype_mat(snp_ind:snp_ind+chr_num_snps_vec(j)-1, cur_sample) = DispStruct.Chrom{j}.Genotypes;
            average_copy_num_mat(snp_ind:snp_ind+chr_num_snps_vec(j)-1, cur_sample) = DispStruct.Chrom{j}.Data(:,3);
            copy_num_mat(snp_ind:snp_ind+chr_num_snps_vec(j)-1, cur_sample) = DispStruct.Chrom{j}.Data(:,1)+...
                DispStruct.Chrom{j}.Data(:,2);
            snp_ind = snp_ind + chr_num_snps_vec(j);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HMMOutStruct.ModelParams = HMM_MODEL; % Save also the Model we've learned, with very long PLACE_M
    if(~exist(hmm_out_dir)) % Make the directory if needed
        eval(['mkdir ' hmm_out_dir]);
    end

    % Save only if score is improved ..
    if(exist(fullfile(user_dir, hmm_out_file_name)))
        OLD = load(fullfile(user_dir, hmm_out_file_name));
        if(HMMOutStruct.LogScore >= OLD.HMMOutStruct.LogScore)
            save(fullfile(user_dir, hmm_out_file_name), 'HMMOutStruct');
        end
        clear OLD;
    else
        save(fullfile(user_dir, hmm_out_file_name), 'HMMOutStruct');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set the output file names:
    OutputFilesNames.disp_files{cur_sample+1} = [display_file_name '.mat'];
    OutputFilesNames.hmm_out_files{cur_sample+1} = [hmm_out_file_name '.mat'];

    %%    close ([1:5]); % close figures but not the GUI window ...

    ttt_MIX = cputime - ttt_MIX
    ttt = cputime - ttt
    fprintf(fid,'sample %ld  rest time %lf\n',cur_sample, ttt);


end % loop on samples
%remove blank SNPs
idx=find(strcmp(deblank(data_snp_ids),''));
average_copy_num_mat(idx,:)=[];
copy_num_mat(idx,:)=[];
chr_num_snps_vec = chr_num_snps_vec - num_of_blank_SNPs;
genotype_mat(idx,:)=[];
data_snp_ids(idx)=[];
%%%%%%%%%%%%%%%%%
ret_file = fullfile(user_dir, 'display', ['all_samples_hmm_out_' chip_type '.mat']);
sample_names = SampleNames;
save(ret_file, 'average_copy_num_mat',  'copy_num_mat', 'chr_num_snps_vec', 'genotype_mat', 'data_snp_ids', 'sample_names');
OutputFilesNames.disp_files{1} = fullfile('display', 'AllSamplesAverage_disp.mat'); % Add also files names fore averages
OutputFilesNames.hmm_out_files{1} = fullfile('hmm_out', 'AllSamplesAverage.mat');

% Calculate the average - we want this to create also a disp file !!!
if(num_samples > 1)
    AverageOverAllSamples(user_dir, SampleNames, chip_type);
end

fclose(fid);

% save genotyping info in normalizarion and hmm summary file
add_genotype_info_to_norm_and_hmm_summary_file(user_dir, chip_type, genotype_mat, sample_names);


% Compare them with Affymetrix's genotypes:
% Syntax of calling the comparison function: 
%% HapMapCompareAffyAndHMMGenotypes(SampleNames, chip_type, user_dir, ALL_MAT, ...
%%    comp_hap_affy, comp_hap_hmm, comp_affy_hmm, use_strand_vec, hapmap_dir, hapmap_population, hapmap_version)

comp_hap_affy = 1; comp_hap_hmm = 1; comp_affy_hmm = 1; use_strand_vec = [0 1 1];
%% TTT_error = cputime;
% ErrorStructWithNoCorrsOld = HapMapCompareAffyAndHMMGenotypes(SampleNames, chip_type, [user_dir '\'], ...
%    comp_hap_affy, comp_hap_hmm, comp_affy_hmm, use_strand_vec);
%% use_strand_vec = [0 0  ];
%% ErrorStructWithNoCorrsLeukemia = HapMapCompareAffyAndHMMGenotypes(SampleNames, chip_type, [user_dir '\'], ...
%% comp_hap_affy, comp_hap_hmm, comp_affy_hmm, use_strand_vec);
%% FirstSample{1} = SampleNames{1};
%% ErrorStructLeukemia = HapMapCompareAffyAndHMMGenotypes(SampleNames, chip_type, [user_dir '\'], 0, 0, HMMParamsStruct.use_affy_genotypes, [0 0 1]);
%% save 'ErrorStructWithCorrs.mat' 'ErrorStructWithCorrs';
%% cputime - TTT_error

ttt_all = cputime - ttt_all

TTT = 999;