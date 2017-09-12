% A function for learning the HMM models from the data. Here is what the
% function does:
%
% Need to fill !!!
%
function [HMM_MODEL_SAVED, Viterbi_Path, Gamma_Probs, HMM_CHROM_KL_DIST, PatChromMeanMatrix, PatChromStdMatrix] = ...
    LearnHMMChromFromData(full_path_data_file_name, full_path_model_and_output_file_name, DataType, ...
    ChromosomesToRun, sample_name, HMM_x_dim, HMM_y_dim, ...
    learn_model_EM, find_path_Viterbi, compute_chromosomes_distance_matrix, use_locations, ...
    do_center_norm, do_fold_change, do_determine_x_dim, do_median_flag, ...
    LD, snp_ids, rs_ids, chr_loc_vec, chr_vec, strand, chip_type, user_dir, genome_assembly, ...
    num_EM_iters, num_EM_starting_points, EM_tolerance, do_load)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HMMOutStruct = [];
hmm_out_file_name = [user_dir  '/hmm_out/' sample_name '_' chip_type '_hmm_out']; 
RelevantInds = [1:6 17:22 33:38 49:54 65:70 81:86]; % This is good for the case 
for i=1:24
    HMMOutStruct.SNPsIDs{i} = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Various data types
MRNA_EXP = 0; SNP_CHIPS=1;
DISCRETE = 0; GAUSSIAN = 1;

% load data from some chromosome
path(path,'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom');  % weizmann
path(path,'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\src');  % weizmann
path(path,'C:\Weizmann\HMM_ido_2004_02_16\hmm_chromosome\hmm_chrom\hmm_chrom');  % home
path(path,'C:\Weizmann\HMM_ido_2004_02_16\hmm_chromosome\hmm_chrom');  % home
path(path,'C:\Weizmann\hmm_chrom\hmm_chrom');  % home

% Load the needed data of the chip
load(full_path_data_file_name);

% In the new version, we need to do some conversions of the data ...
if(DataType == SNP_CHIPS) % Work for now only on the current SNP chip
    SPECIAL_MODELS_FLAG = 1; % Must have the special SNPs model
    
    snp_id_str = ['snp_id_' lower(chip_type)];
    eval(['snp_id_chip = ' snp_id_str ';']);
    [SnpsNames I J] = intersect(snp_ids, snp_id_chip);
    rs_ids = rs_ids(I); snp_ids2 = snp_ids(I); % Keep also the snp ids. 
    HMM_genes = snp_id_chip;
    HMM_locations = chr_loc_vec(I); % Currently we don't know the locations
    HMM_samples = {}; HMM_samples{1} = sample_name; % Pick one of the samples
    HMM_ref_labels = 1;
    HMM_chromosome_arm = chr_vec(I); % Problem ! here we've got only the chromosome and not the arm!
    sample_ratio_str = ['allele_ratio_vec_' lower(chip_type)];
    sample_copy_num_str = ['copy_num_vec_' lower(chip_type)];
    sample_genotype_call_str = ['genotype_vec_' lower(chip_type)];
    eval([sample_ratio_str '= min(' sample_ratio_str ', 9999999999);']);
    eval(['HMM_data = ' sample_copy_num_str './ (' sample_ratio_str ' + 1);']);
    eval(['HMM_data2 = HMM_data .* ' sample_ratio_str ';']);
    HMM_data = HMM_data(J); HMM_data2 = HMM_data2(J);
    TOL = 0.00000000000001;
    eval(['max_error = max( (HMM_data+HMM_data2-' sample_copy_num_str ').^2)']);
end


% Here's what we got : (old version)
% HMM_chromosome_arm   - which arm (e.g. 7q) every gene lies on
% HMM_data             - samples (expression data matrix ) of all genes
% HMM_genes            - labels of all genes
% HMM_locations        - locations (in nucleotides from chromosome starts of each gene on chromosome
% HMM_samples          - labels of all samples
% HMM_ref_labels       - labels saying if this sample is in the reference ('Normal') group or in the set we want to check ('Cancer')

plot_vec = 'byrgmck';  % Colors for plotting
% Vector denoting the start of the q-arm in each chromosome
first_on_q =[142604331, 95203899, 94912966, 52697919, 50705877, 62387336, 63850117, 48697951,  66551451, 41965097, 55478369, ...
    39709262, 17055261, 18814688, 18962068, 46771070, 25766888, 17462263, 34390067, 30817204, 14665357, 14608115];

epsilon = 0.000000001;
num_chroms = 22;
do_log = 0; % Flag saying if to perform log transform
% do_center_norm = 1; % Flag saying if to do centering and normalization
place_flag = DataType; % give a different distribution for every place. These are based on chromosome 3
% which is divided to 10 patients with LOH and 10 patients with NLOH
cv = 0;   % cross-validation technique - this is the number of samples we learn on each time.
% 0 means no cv - training and testing on same data (which is wrong!!)

num_genes = length(HMM_genes); num_samples = length(HMM_samples);

% First do fold-change if neccessary. We will keep also the healthy ones !!!!!!
if(do_fold_change)
    HMM_take_indexes = find(HMM_ref_labels);
    HMM_ref_indexes = find(HMM_ref_labels==0);
    if(do_median_flag == 0)
        HMM_ref_mean = mean(HMM_data(:,HMM_ref_indexes), 2);  % Must be first !!!
    else
        HMM_ref_mean = median(HMM_data(:,HMM_ref_indexes), 2);
    end
    num_samples = length(HMM_samples); % Correct the number of samples
    % Now make the fold-change !!! and take the log !!
    HMM_data = HMM_data ./ repmat(HMM_ref_mean, 1, num_samples);
    HMM_data = log(HMM_data);
    % Make a plot of the means of patients log-fold-change at each chromosome
    PatChromMeanMatrix = zeros(num_chroms, num_samples);
    PatChromStdMatrix = zeros(num_chroms, num_samples);
end

% Get rid of this stupid cell format, and get data as numbers
if(DataType == SNP_CHIPS)
    HMM_chromosome = HMM_chromosome_arm; % Here we already have the chromosomes only
% %     HMM_chromosome = (char(HMM_chromosome_arm)); % Here we already have the chromosomes only
% %     XXX = find(HMM_chromosome == 'X');
% %     HMM_chromosome(XXX,1) = '2'; % Encode the X chromosome
% %     HMM_chromosome(XXX,2) = '3'; % Encode the X chromosome
% %     III = intersect( find(HMM_chromosome(:,1) == ' '), find(HMM_chromosome(:,2) == ' ') );
% %     HMM_chromosome(III,1) = '0'; % Encode the chromosomes not known at all
% %     HMM_chromosome = str2num(HMM_chromosome);
else
    HMM_chromosome = zeros(1,num_genes);
    for i=1:num_genes
        HMM_chromosome(i) = str2num(HMM_chromosome_arm{i}(find((HMM_chromosome_arm{i} <= '9')&(HMM_chromosome_arm{i} >= '0')))); % get only the numerical digits
    end
end
num_genes_in_chrom = zeros(1,num_chroms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now do loading of the models
% When and what we load and what do we compute every time???
cur_str = 'hmm_models';
if(use_locations)
    cur_str = [cur_str '_use_loc'];
end

if(cv > 0)  % here cross-validation
    cur_str = [cur_str '_cv_out' num2str(cv_test_set(1)) '-' num2str(cv_test_set(end))];
end

do_load_is = do_load

if(do_load)
    load(full_path_model_and_output_file_name); % Note : This also contains the HMM_MODEL_Details variable which is now problematic !
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HMM_MODEL = {};

HMM_chrom_loc_mat_copy = {};
HMM_chrom_data_mat_copy = {};   % Now we set the data also in cell file !!
HMM_chrom_data_mat_copy2 = {};

all_data_mat = []; all_data_A_mat = []; all_data_B_mat = []; all_loc_mat = [];

% Here do the Viterbi and Forward/Backward part
for i=ChromosomesToRun %%% num_chroms % go over all chromosomes

    NowDoDataCollectingChromosome = i
    
    % find which genes are on it
    indexes_vec = find(HMM_chromosome == i);

    [JointSNPs III JJJ] = intersect(LD{i}.RS_IDs, rs_ids(indexes_vec)); % need to do further intersection since some SNPs are ? have the same rs_ids ?

    HMMOutStruct.SNPsIDs{i} = snp_ids2(indexes_vec); 
    indexes_vec = indexes_vec(JJJ);
    MissInds = setdiff( [1:length(LD{i}.RS_IDs)], III );

    % Take the sub-matrix of genes on this chromsome
    HMM_chrom_data_mat_copy{i} = zeros(length(LD{i}.RS_IDs),1);
    HMM_chrom_data_mat_copy2{i} = zeros(length(LD{i}.RS_IDs),1);
    HMM_chrom_loc_mat_copy{i} = zeros(length(LD{i}.RS_IDs),1);
    HMM_chrom_data_mat_copy{i}(III) = HMM_data(indexes_vec,:);
    HMM_chrom_data_mat_copy2{i}(III) = HMM_data2(indexes_vec,:);
    HMM_chrom_loc_mat_copy{i}(III) = HMM_locations(indexes_vec);
    for j=MissInds
        HMM_chrom_data_mat_copy{i}(j) = ...
            (HMM_chrom_data_mat_copy{i}(max(j-1,1)) + HMM_chrom_data_mat_copy{i}(min(j+1,length(LD{i}.RS_IDs))))/2;
        HMM_chrom_data_mat_copy2{i}(j) = ...
            (HMM_chrom_data_mat_copy2{i}(max(j-1,1)) + HMM_chrom_data_mat_copy2{i}(min(j+1,length(LD{i}.RS_IDs))))/2;
        HMM_chrom_loc_mat_copy{i}(j) = ...
            (HMM_chrom_loc_mat_copy{i}(max(j-1,1)) + HMM_chrom_loc_mat_copy{i}(min(j+1,length(LD{i}.RS_IDs))))/2;
    end

    all_data_mat = [all_data_mat HMM_chrom_data_mat{i}'+HMM_chrom_data_mat2{i}'];
    all_data_A_mat = [all_data_A_mat HMM_chrom_data_mat{i}'];
    all_data_B_mat = [all_data_B_mat HMM_chrom_data_mat2{i}'];
    all_loc_mat = [all_loc_mat HMM_chrom_loc_mat_copy{i}'];

end

sm_all_data_mat = smooth(all_data_mat, 100);
sm_all_data_A_mat = smooth(all_data_A_mat, 30);
sm_all_data_B_mat = smooth(all_data_B_mat, 30);
% First run a 'primitive' MOG just to see if everyhthing fits
% ttt = cputime; [P_s,M_s,S_s, Dim_s, LogLike_s]=MixtureOfGaussiansFindModelDimension(sm_all_data_mat,5,50, 10);  cputime - ttt
% % Now compensate for the other levels not seen here ...
% [P_max P_max_ind] = max(P_s);
% SS = zeros(1,5); SS(3) = S_s(P_max);
% for i=1:Dim_s
%     SS(i-P_max_ind+3) = S_s(i);
% end


% Learn a model for all the chromosomes of this patient
% here learn the model parameters
if(learn_model_EM)

    % Dummy variables for calling c-function
    mean_vec_rep = zeros(num_genes_in_chrom(i)*num_samples,HMM_x_dim);
    std_vec_rep = zeros(num_genes_in_chrom(i)*num_samples,HMM_x_dim);

    % Note : Here we try to learn the SIZE of the model (the number of
    % states X can attain !!
    if(do_determine_x_dim)
        max_x_dim = HMM_x_dim; % Now this is the maximal value, while the actual value might be smaller
        [K_min, pen_dist_vec, params_cell, logscores] = ...
            HMM_K_by_penalized_dist_EM(HMM_chrom_data, [1:max_x_dim], num_EM_iters, num_EM_starting_points, EM_tolerance);

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

        %  HMM_x_dim = K_min;
    else  % Here we know the model size
        K_min = HMM_x_dim; % HMM_x_dim
        Ksize = HMM_x_dim;
    end  % determine x_dim

    % Note : We assume a stationary distribution here so PI is
    % determined accordingly
    % Learn Again ! (This should be unneccessary) - we learn a different
    % model for each sample, but the model is the SAME for all the chromosomes !!!! 

    [HMM_MODEL.PI HMM_MODEL.M HMM_MODEL.N HMM_MODEL.MEW ...
        HMM_MODEL.SIGMA HMM_MODEL.LogScore] = ...
        TrainHMMFromDataEMMatlab(all_data_A_mat, all_data_B_mat, HMM_chrom_loc, use_locations, K_min, HMM_y_dim, ...
        0, do_fold_change, mean_vec_rep, std_vec_rep, ...
        SPECIAL_MODELS_FLAG, eye(Ksize), 0, ... % Here we put the upper-bounds. Currently they're not used !!!
        num_EM_iters, num_EM_starting_points, EM_tolerance);

    %
    %             /* Parameters should be :
    % 	1. Expression vector - 1 * (num_samples * num_genes)
    % 	2. Chromosomal location vector  - same
    % 	3. Missing data flag (says if to use locations ... recommended zero)
    % 	4. X dimension required   (number of hidden states)
    % 	5. Y dimension (mixtures) required   (number of mixtures of gaussians)
    % 	6. place flag (if to use different distribution for each place) (recommended zero)
    % 	7. fold change flag (if to bound the mews by one)  (recommended zero)
    % 	8. mean vec (currently it is given and not learned)  (dummy, unless place flag is one)
    % 	9. std vec (currently it is given and not learned)  (dummy, unless place flag is one)
    % 	10. number of EM iterations   (50-100)
    % 	11. Number of EM starting points   (5-10)
    % 	12. EM tolerance     (10^(-6))
    % 	*/


    % Function returns :
    %  PI - stationary dist. of hidden variable
    %  M - Transition mat.
    %  N - Mixture matrix
    %  MEW - mean matrix
    %  SIGMA - std. dev. matrix
    %  LogScore - Log Likelihood of data given the returned parameters


    CurrentLogScore = HMM_MODEL.LogScore
    if(do_load)
        copy_flag = 0;
        if(copy_flag == 1)
            HMM_MODEL_SAVED = HMM_MODEL;
            HMM_chrom_loc_mat = HMM_chrom_loc_mat_copy;
            HMM_chrom_data_mat = HMM_chrom_data_mat_copy;
            HMM_chrom_data_mat2 = HMM_chrom_data_mat_copy2;
        end
    else
        HMM_MODEL_SAVED = HMM_MODEL;  % Here we save everything because this is the first time ...
        copy_flag = 1;
    end
else % Do not learn parameters
    K_min = HMM_x_dim; % HMM_x_dim
    Ksize = HMM_x_dim;
    copy_flag = 1;
end   % finished learning part



% Here do the Viterbi and Forward/Backward part
for i=ChromosomesToRun %%% num_chroms % go over all chromosomes

    NowDoViterbiChromosome = i

    num_genes_in_chrom(i) = length(LD{i}.RS_IDs);
    % New ! Sort the genes such that locations are in an increasing order !
    [ HMM_chrom_loc_mat_copy{i} perm ] = sort(HMM_chrom_loc_mat_copy{i});
    HMM_chrom_data_mat_copy{i} = HMM_chrom_data_mat_copy{i}(perm,:);
    HMM_chrom_data_mat_copy2{i} = HMM_chrom_data_mat_copy2{i}(perm,:);
    % Reshape to a vector form
    HMM_chrom_data = reshape(HMM_chrom_data_mat_copy{i}, 1, num_genes_in_chrom(i)*num_samples);
    HMM_chrom_data2 = reshape(HMM_chrom_data_mat_copy2{i}, 1, num_genes_in_chrom(i)*num_samples);
    HMM_chrom_loc = (repmat(HMM_chrom_loc_mat_copy{i}, num_samples, 1))';


    % What do we copy here??
    if(copy_flag == 1)
        if(DataType == SNP_CHIPS) % Fill all the needed details (without learning ... )
            derich = 1.0/120.0; % Relative Derichlet correction
            % Determine parameters ...
            P_chrom_copy_number_change = 0.001;  % Probability to change the copy number
            HMM_MODEL_SAVED{i}{K_min}.SIGMA = [0.0778    0.2868    0.4103    0.3611    0.6931]; % From MOG fitting
            %%%%%%  0.645 * ones(1,2*HMM_x_dim-1); % Increase sigma ... remember that we have more than x_dim
            HMM_MODEL_SAVED{i}{K_min}.N = ones(2*HMM_x_dim-1,1);
            HMM_MODEL_SAVED{i}{K_min}.PI = ones(HMM_x_dim,1) ./ HMM_x_dim;
            HMM_MODEL_SAVED{i}{K_min}.M = eye(HMM_x_dim) .* (1-2*P_chrom_copy_number_change);
            for j=1:HMM_x_dim-1
                HMM_MODEL_SAVED{i}{K_min}.M(j+1,j) = P_chrom_copy_number_change;
                HMM_MODEL_SAVED{i}{K_min}.M(j,j+1) = P_chrom_copy_number_change;
            end
            HMM_MODEL_SAVED{i}{K_min}.M(1,1) = 1-P_chrom_copy_number_change;
            HMM_MODEL_SAVED{i}{K_min}.M(HMM_x_dim,HMM_x_dim) = 1-P_chrom_copy_number_change;
            HMM_MODEL_SAVED{i}{K_min}.MEW = [ 0.1422 0.7262 1.5627 2.3612 2.4599]'; % From MOG fitting %%%%[0.2 0.8 1.8 2.45 3]';  % Note: Mew here is of size 2*dim-1 !!!
            HMM_MODEL_SAVED{i}{K_min}.PLACE_FLAG=1;
            HMM_MODEL_SAVED{i}{K_min}.SPECIAL_MODELS_FLAG=1;
            % Calc the LD transition matrices - Transfer from joint to conditional probabilities:
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M = zeros(4,length(LD{i}.PairMat(:,1)));
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M(1,:) = (LD{i}.PairMat(:,1)+derich) ./ (LD{i}.PairMat(:,1)+LD{i}.PairMat(:,2)+2*derich);
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M(3,:) = 1-HMM_MODEL_SAVED{i}{K_min}.PLACE_M(1,:);  % Switch 2 and 3 for the C function
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M(2,:) = (LD{i}.PairMat(:,3)+derich) ./ (LD{i}.PairMat(:,3)+LD{i}.PairMat(:,4)+2*derich);
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M(4,:) = 1-HMM_MODEL_SAVED{i}{K_min}.PLACE_M(2,:);

            % Calc the 'singleton' LD transition matrices. Take into account only frequencies, and not the pairwise LD correlations
            HT = strand(I);
            HT_SIGNS = zeros(1,length(HT));
            HT_SIGNS(strmatch('+', HT)) = 1; % '+' strand is 1
            CHR_HT_SIGNS = zeros(1,length(LD{i}.FreqVec));
            CHR_HT_SIGNS(III) = HT_SIGNS(indexes_vec);
            %%% LD_corrected_FreqVec = 2 .* LD{i}.FreqVec .* CHR_HT_SIGNS' - LD{4}.FreqVec - CHR_HT_SIGNS' + 1;

            HMM_MODEL_SAVED{i}{K_min}.PLACE_M2 = zeros(4,length(LD{i}.PairMat(:,1)));
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M2(1,:) = (LD{i}.FreqVec(2:end)+derich) ./ (1+2*derich);
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M2(3,:) = 1-HMM_MODEL_SAVED{i}{K_min}.PLACE_M2(1,:);  % Switch 2 and 3 for the C function
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M2(2,:) = (LD{i}.FreqVec(2:end)+derich) ./ (1+2*derich);;
            HMM_MODEL_SAVED{i}{K_min}.PLACE_M2(4,:) = 1-HMM_MODEL_SAVED{i}{K_min}.PLACE_M2(2,:);

            % Don't use all this LD data at the beginning
            %%%%HMM_MODEL_SAVED{i}{K_min}.PLACE_M = HMM_MODEL_SAVED{i}{K_min}.PLACE_M*0 + 0.5;

            % Problem: We may have some (very few) SNPs which are not in the HAPMAP but on the chip.
            % We want to insert them also so we need to 'strech' the M's (or 'remove' the data).

            % Copy also the data in the SNP case
            HMM_chrom_loc_mat{i} = HMM_chrom_loc_mat_copy{i};
            HMM_chrom_data_mat{i} = HMM_chrom_data_mat_copy{i};
            HMM_chrom_data_mat2{i} = HMM_chrom_data_mat_copy2{i};
        end
        % This actualy does the 'classification'
        if(find_path_Viterbi)
            Viterbi_Path_copy{i} = zeros(num_samples, num_genes_in_chrom(i) );  % how many genes are in the current chromosome
            %%%            Gamma_Probs_copy{i} = zeros(num_samples, HMM_MODEL_Details.KMin(i), num_genes_in_chrom(i) );  % how many genes are in the current chromosome
            Gamma_Probs_copy{i} = zeros(num_samples, 256, num_genes_in_chrom(i) );  % how many genes are in the current chromosome

            % Dummy again ...
            mean_vec = zeros(1,num_genes_in_chrom(i));
            std_vec = 1+zeros(1,num_genes_in_chrom(i));
            for j = 1:num_samples
                HMM_MODEL_SAVED{i}{K_min}.Y_TYPE = GAUSSIAN;
                HMM_MODEL_SAVED{i}{K_min}.y_dim = 1; % No mixture, only one gaussian ...
                % See if the model is
                %%%   [VVV GGG] = ...
                p = 1;
                HMM_MODEL_SAVED{i}{K_min}.M2 = p .* HMM_MODEL_SAVED{i}{K_min}.M + (1-p) .* eye(HMM_x_dim);
                %                 HMM_MODEL_SAVED{i}{K_min}.SIGMA = 0*HMM_MODEL_SAVED{i}{K_min}.SIGMA+0.645;
                start_ind = 1; end_ind = length(HMM_chrom_data_mat_copy{i});
                Y = HMM_chrom_data_mat_copy{i}(start_ind:end_ind,j);
                Y2 = HMM_chrom_data_mat_copy2{i}(start_ind:end_ind,j);
                tmp_Y = Y; Y = Y .* CHR_HT_SIGNS(start_ind:end_ind)' + Y2 .* (1-CHR_HT_SIGNS(start_ind:end_ind)');
                Y2 = tmp_Y .* (1-CHR_HT_SIGNS(start_ind:end_ind)') + Y2 .* CHR_HT_SIGNS(start_ind:end_ind)'; % Switch Y and Y2 where strand is negative
                L = HMM_chrom_loc_mat_copy{i}(start_ind:end_ind);

                [Viterbi_Path_copy{i}(j,:) Gamma_Probs_copy{i}(j,:,:)] = ...
                    FindBestPathViterbi(Y, Y2, L, ...
                    use_locations, HMM_MODEL_SAVED{i}{K_min}.PI, HMM_MODEL_SAVED{i}{K_min}.M2, HMM_MODEL_SAVED{i}{K_min}.N, ...
                    HMM_MODEL_SAVED{i}{K_min}.MEW', HMM_MODEL_SAVED{i}{K_min}.SIGMA, HMM_MODEL_SAVED{i}{K_min}.PLACE_FLAG, ...
                    HMM_MODEL_SAVED{i}{K_min}.SPECIAL_MODELS_FLAG, mean_vec, std_vec, HMM_MODEL_SAVED{i}{K_min}.PLACE_M);  % add the place variables all_chr_exp_arr(j,:);
                do_couples = 1;
%%                [V_VEC_alpha_genotype V_VEC_beta_genotype V_VEC_alpha_copynumber{i} V_VEC_beta_copynumber{i} ...
%%                    V_VEC_copy V_VEC_A V_VEC_B] = ...
                 [VitStruct ProbsStruct] = ...    
                    GetBestMarginalPredictions(reshape(Gamma_Probs_copy{i}(j,:,:), 256, length(Gamma_Probs_copy{i}(j,:,:))), HMM_MODEL_SAVED{i}{K_min}, do_couples);

                Viterbi_Path_copy{i} = Viterbi_Path_copy{i}'; %%Gamma_Probs{i} = Gamma_Probs{i}';
            end % loop on patients
            Gamma_Probs{i} = Gamma_Probs_copy{i};
            Viterbi_Path{i} = Viterbi_Path_copy{i};
        end  % Viterbi
    end % copy flag
end  % loop on chromosomes (i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Transfer V_VEC to its 'real' values
% [V_VEC_alpha_genotype V_VEC_beta_genotype V_VEC_alpha_copynumber V_VEC_beta_copynumber ...
%     V_VEC_copy V_VEC_A V_VEC_B] = ...
%     GetBestMarginalPredictions(reshape(Gamma_Probs{i}(1,:,:), 256, size(Gamma_Probs{i}, 3)), K_min);  % Need to change whole input here ...


% % % % % % % % % V_VEC_alpha_genotype = bitget(Viterbi_Path{i},1);
% % % % % % % % % V_VEC_beta_genotype = bitget(Viterbi_Path{i},5);
% % % % % % % % % V_VEC_alpha_copynumber = bitget(Viterbi_Path{i},2) + 2*bitget(Viterbi_Path{i},3) + 4*bitget(Viterbi_Path{i},4);
% % % % % % % % % V_VEC_beta_copynumber = bitget(Viterbi_Path{i},6) + 2*bitget(Viterbi_Path{i},7) + 4*bitget(Viterbi_Path{i},8);
% % % % % % % % % V_VEC_copy = V_VEC_alpha_copynumber+V_VEC_beta_copynumber;
% % % % % % % % % V_VEC_A = V_VEC_alpha_genotype.*V_VEC_alpha_copynumber + ...
% % % % % % % % %     V_VEC_beta_genotype.*V_VEC_beta_copynumber;
% % % % % % % % % V_VEC_B = (1-V_VEC_alpha_genotype).*V_VEC_alpha_copynumber + ...
% % % % % % % % %     (1-V_VEC_beta_genotype).*V_VEC_beta_copynumber;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
[LLL LLL_I LLL_J] = intersect(HMM_chrom_loc_mat_copy{i}, HMM_locations); % Get the locations of the SNPs we still have
eval(['[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(' sample_genotype_call_str '(LLL_J));']);
%%%%[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(HD78_9_diag_calls_hind(LLL_J)); % Get the Genotypes
figure; hold on;
plot(HMM_chrom_data_mat_copy{i}(LLL_I(AA_ind1),j),HMM_chrom_data_mat_copy2{i}(LLL_I(AA_ind1),j), '.');
plot(HMM_chrom_data_mat_copy{i}(LLL_I(AB_ind1),j),HMM_chrom_data_mat_copy2{i}(LLL_I(AB_ind1),j), 'r.');
plot(HMM_chrom_data_mat_copy{i}(LLL_I(BB_ind1),j),HMM_chrom_data_mat_copy2{i}(LLL_I(BB_ind1),j), 'g.');
plot(HMM_chrom_data_mat_copy{i}(LLL_I(no_call_ind1),j),HMM_chrom_data_mat_copy2{i}(LLL_I(no_call_ind1),j), 'm.');
legend('AA', 'AB', 'BB', 'no call');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


% In the first time we copy everything and save !!!!
if(do_load == 0)
    HMM_chrom_loc_mat = HMM_chrom_loc_mat_copy;
    HMM_chrom_data_mat = HMM_chrom_data_mat_copy;
    HMM_chrom_data_mat2 = HMM_chrom_data_mat_copy2;
    Gamma_Probs = Gamma_Probs_copy;
    Viterbi_Path = Viterbi_Path_copy;
end

% Saving !!!
saving = full_path_model_and_output_file_name

% New !! add some parameters to the file !!!
save(full_path_model_and_output_file_name, 'HMM_MODEL_SAVED', 'Viterbi_Path', 'Gamma_Probs', ...
    'HMM_CHROM_KL_DIST', 'num_genes_in_chrom', 'HMM_chrom_loc_mat', 'HMM_chrom_data_mat',  'HMM_chrom_data_mat2', ...
    'num_genes_in_chrom', 'PatChromMeanMatrix', 'PatChromStdMatrix');


% Newest : Save in the format required by: 1. The displaying program
%                                          2. Further multi-SNP analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. The displaying program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Upperbound = 2.5; % amplifications
Lowerbound = 1.5; % deletions
display_file_name = [user_dir  '/display/' sample_name '_' chip_type '_disp']; 
DispStruct = [];
DispStruct.SampleName = sample_name;
DispStruct.ChipName = chip_type;
DispStruct.GenomeBuild = genome_assembly;
for i=1:24
    DispStruct.Chrom{i} = [];
end
for i=ChromosomesToRun
    DispStruct.Chrom{i}.Data = [V_VEC_alpha_copynumber{i} V_VEC_beta_copynumber{i}];
    DispStruct.Chrom{i}.Locs = HMM_chrom_loc_mat{i};
%%    DispStruct.Chrom{i}.Segments = ??;%%% (So we do here also the segment analysis ??)
    AmpInds = (V_VEC_alpha_copynumber{i}+V_VEC_beta_copynumber{i} > Upperbound)'; % Find the segments positions
    [AmpIndsStarts AmpIndsEnds] = GetSegmentsFromIndexes(AmpInds);
    DelInds = (V_VEC_alpha_copynumber{i}+V_VEC_beta_copynumber{i} < Lowerbound)';
    [DelIndsStarts DelIndsEnds] = GetSegmentsFromIndexes(DelInds);
    AllIndsStarts = [AmpIndsStarts DelIndsStarts];
    AllIndsEnds = [AmpIndsEnds DelIndsEnds];
    [AllIndsStarts SortPerm] = sort(AllIndsStarts); AllIndsEnds = AllIndsEnds(SortPerm); % Sort segments
    DispStruct.Chrom{i}.Segments(:,1) = HMM_chrom_loc_mat{i}(AllIndsStarts);
    DispStruct.Chrom{i}.Segments(:,2) = HMM_chrom_loc_mat{i}(AllIndsEnds);
    V_VEC_cumsum = [0 cumsum(V_VEC_alpha_copynumber{i}+V_VEC_beta_copynumber{i})'];
    DispStruct.Chrom{i}.Segments(:,3) = (V_VEC_cumsum(AllIndsEnds+1) - V_VEC_cumsum(AllIndsStarts)) ./ ...
        (AllIndsEnds-AllIndsStarts+1); % Find the average of each semgment
end
save display_file_name DispStruct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2. Further multi-SNP analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hmm_out_file_name = [user_dir  '/hmm_out/' sample_name '_' chip_type '_hmm_out']; 
RelevantInds = [1:6 17:22 33:38 49:54 65:70 81:86]; % This is good for the case 
for i=1:24
    HMMOutStruct.SNPsProbsMat{i} = [];
end
for i=ChromosomesToRun
    HMMOutStruct.SNPsProbsMat{i} = reshape(Gamma_Probs{i}(1,RelevantInds,:), 4*HMM_x_dim^2, size(Gamma_Probs{i},3)); 
end
save hmm_out_file_name HMMOutStruct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
