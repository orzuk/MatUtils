% Written by Or Zuk 5/2007
%
% Find what happens if we generate data from a given model and try to learn
% it again. This function gives benchmarks for recovering the true
% sequence, the true model parameters etc. - Everything is caclulated
% when we know exactly (and control) the process which generated the data
function [COMP_HMM_MODELS] = DetermineModelStability( HMM_MODEL, num_points_vec, num_iters, ...
    num_EM_iters, num_EM_starting_points, learn_flag, learn_equ_flag)

% Tables for indexing the A and B copy 
B_copy_tab = [ 0, 0, 0, 1, 0, 2, 0, 0, 0, ...
               1, 0, 2, 0, 0, 0, 1, 0, 2, ...
               1, 1, 1, 2, 1, 3, 0, 0, 0, ...
               1, 0, 2, 2, 2, 2, 3, 2, 4];
A_copy_tab = [ 0, 0, 1, 0, 2, 0, 0, 0, 1, ...
               0, 2, 0, 1, 1, 2, 1, 3, 1, ...
               0, 0, 1, 0, 2, 0, 2, 2, 3, ...
               2, 4, 2, 0, 0, 1, 0, 2, 0]; 
%% Total_copy_tab = A_copy_tab + B_copy_tab;
%% AB_copy_tab = A_copy_tab + 5*B_copy_tab; % here put A in the lsb (base 5) and B in the  msb (another 5)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some default values for parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_couples = 0; %%% 1; % We need to take the coupling into account
joint_flag = 1; %%% 1; % We compute marginal over the joint distribution of both chromosomes
% Parameters for KL distance
KL_dist_seq_len = 2500; % Length of sequence to simulate for caclulating KL
num_KL_dist_iters = 10; % Num. of iterations when generating sequences for calculating KL
MEW_GAP  = max(HMM_MODEL.MEW) - min(HMM_MODEL.MEW); % The gap between the mew's for 'normalizing' the results


% Do many simulations, and each time determine the best learned model
COMP_HMM_MODELS = {}; K = length(num_points_vec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These structures contain the errors of various kinds, which we later plot to check our performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS_ERROR = InitParamsError(K, HMM_MODEL);
CLASS_ERROR = InitClassError(K);
CLASS_ERROR_SINGLETONS = InitClassError(K);
CLASS_ERROR_IID = InitClassError(K);

% Prepare alternative (less informative) LD matrices
HMM_MODEL.PLACE_M_IID = 0*HMM_MODEL.PLACE_M + 0.5;
CurVec = [0.5 0.5]; % We don't know the frequencies for the first vector (do we? look in simulator!!)
HMM_MODEL.PLACE_M_SINGLETONS = HMM_MODEL.PLACE_M;
HMM_MODEL.PLACE_M_SINGLETONS(:,1) = 0.5;
MMM = reshape(HMM_MODEL.PLACE_M, 2, 2, length(HMM_MODEL.PLACE_M));
for i=1:length(HMM_MODEL.PLACE_M)-1
    CurVec = CurVec * MMM(:,:,i);
    HMM_MODEL.PLACE_M_SINGLETONS(1,i+1) = CurVec(1);
    HMM_MODEL.PLACE_M_SINGLETONS(2,i+1) = CurVec(1);
    HMM_MODEL.PLACE_M_SINGLETONS(3,i+1) = CurVec(2);
    HMM_MODEL.PLACE_M_SINGLETONS(4,i+1) = CurVec(2);
end
%% HMM_MODEL.PLACE_M_SINGLETONS(1,:) = mean(HMM_MODEL.PLACE_M(1:2,:));
%% HMM_MODEL.PLACE_M_SINGLETONS(2,:) = mean(HMM_MODEL.PLACE_M(1:2,:));
%% HMM_MODEL.PLACE_M_SINGLETONS(3,:) = mean(HMM_MODEL.PLACE_M(3:4,:));
%% HMM_MODEL.PLACE_M_SINGLETONS(4,:) = mean(HMM_MODEL.PLACE_M(3:4,:));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Big loop - Each time do:
% 1. Simulate data
% 2. (Learn parameters if learn flag is on)
% 3. Run Viterbi to classify using original/learned model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
points_ind = 1;  % Keep index of current number of points
for num_points = num_points_vec % loop over sequence length
    for i=1:num_iters     % loop over number of iterations
        % Simulate sequence
        SAVE_M = HMM_MODEL.M;
        if( size(HMM_MODEL.PLACE_M,2) == 4) % Transpose if needed
            HMM_MODEL.PLACE_M = HMM_MODEL.PLACE_M';
        end
        %% Check some specific modifications to see their effect on performance
        %% HMM_MODEL.M = eye(3); % Don't allow ANY transitions !!!
        %% HMM_MODEL.PI = [0 1 0]; % Force the level at the beginning to be one!!!
        [X_VEC.all_data OUT_VEC OUT_VEC_B] = SimulateSequenceFromModel(HMM_MODEL.PI, HMM_MODEL.M, HMM_MODEL.N, ...
            HMM_MODEL.MEW, HMM_MODEL.SIGMA, max(HMM_MODEL.PLACE_FLAG,HMM_MODEL.SPECIAL_MODELS_FLAG),  ...
            HMM_MODEL.PLACE_M, HMM_MODEL.SPECIAL_MODELS_FLAG, num_points); % Simulate Sequence (assume special flag is on)
        HMM_MODEL.M = SAVE_M;         % Change them back for the learning part
        
        % Correct PI
        [eig_vecs eig_vals] = eig(HMM_MODEL.M');
        [min_val min_ind] = min(diag(eig_vals-1).^2);
        HMM_MODEL.PI = eig_vecs(:,min_ind) ./ sum(eig_vecs(:,min_ind));


        % 'Open' and Transfer X_VEC to its 'meaningful' values
        if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
            X_VEC = GetMeaningfullCopyNumbersAndGenotypes(X_VEC.all_data);
            equ_PI = zeros(1,5);
            for j=0:2
                for k=0:2
                    equ_PI(j+k+1) = equ_PI(j+k+1) + HMM_MODEL.PI(j+1)*HMM_MODEL.PI(k+1);
                end
            end
            ModelLogLikelihoodRate =  sum(HMM_MODEL.M .* log(HMM_MODEL.M)) * HMM_MODEL.PI  + ...
                [0.5 0.5 0.5 0.5] * ( HMM_MODEL.PLACE_M(:,1) .* log(HMM_MODEL.PLACE_M(:,1)) ) + ...
                (equ_PI *  (-0.5*(1 + log(2*pi) + 2*log(HMM_MODEL.SIGMA))'))
        else
            ModelLogLikelihoodRate =  sum(HMM_MODEL.M .* log(HMM_MODEL.M)) * HMM_MODEL.PI  + ...
                (HMM_MODEL.PI' *  (-0.5*(1 + log(2*pi) + 2*log(HMM_MODEL.SIGMA))'))
        end

        % Set all parameters
        LOC_VEC = 1:length(OUT_VEC);
        use_locations = 0; % Don't use distances between SNPs
        HMM_y_dim = 1; HMM_x_dim = length(HMM_MODEL.PI);
        do_fold_change = 0; % Note : Fold Change was already done !!!!
        mean_vec_rep = zeros(length(OUT_VEC),HMM_x_dim); std_vec_rep = zeros(length(OUT_VEC),HMM_x_dim);
        EM_tolerance = 0.000000001;

        % New: Compute loglikelihood of generated data given the model:
        use_x_flag = 1; % This flag says if we want to calculate log(Pr(X,Y)) or just log(Pr(Y))
        L=1:length(OUT_VEC); % Not important locations
        DataLogLikelihood = ...
            ComputeHMMLogLikelihood(X_VEC.all_data, OUT_VEC, OUT_VEC_B, use_x_flag, L, ...
            use_locations, HMM_MODEL.PI, HMM_MODEL.M, HMM_MODEL.N, ...
            HMM_MODEL.MEW', HMM_MODEL.SIGMA, max(HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.PLACE_FLAG), ...
            HMM_MODEL.SPECIAL_MODELS_FLAG, mean_vec_rep, std_vec_rep, ...
            HMM_MODEL.PLACE_M);

%%        x_loglike = 0
%%        for j=1:99
%%            x_loglike = x_loglike + log(HMM_MODEL.M(X_VEC.alpha_copy(j)+1,X_VEC.alpha_copy(j+1)+1)) + ...
%%                log(HMM_MODEL.M(X_VEC.beta_copy(j)+1,X_VEC.beta_copy(j+1)+1));
%%        end
%%        x_loglike
        
        % Now Compare the log-likelihood we've got to what we'd expect
        DataLogLikelihoodRate = DataLogLikelihood / length(OUT_VEC)

        % Make an 'equivalent' model which holds information only on the copy number
        if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
            equ_x_dim = HMM_MODEL.x_dim*2-1; % Allow to add the A and B copy numbers
            equ_M = zeros(equ_x_dim);  % Build the 'equivalence' M markov matrix
            prod_M = zeros(HMM_MODEL.x_dim.^2);
            for i1=1:HMM_MODEL.x_dim
                for i2=1:HMM_MODEL.x_dim
                    for j1=1:HMM_MODEL.x_dim
                        for j2=1:HMM_MODEL.x_dim
                            prod_M(i1+HMM_MODEL.x_dim*(i2-1),j1+HMM_MODEL.x_dim*(j2-1)) = ...
                                HMM_MODEL.M(i1,j1)* HMM_MODEL.M(i2,j2);
                            equ_M(i1+i2-1,j1+j2-1) = equ_M(i1+i2-1,j1+j2-1) + ...
                                prod_M(i1+HMM_MODEL.x_dim*(i2-1),j1+HMM_MODEL.x_dim*(j2-1));
                        end
                    end
                end
            end
            equ_norm_vec = [1:HMM_MODEL.x_dim HMM_MODEL.x_dim-1:-1:1];
            equ_M = equ_M ./ repmat(equ_norm_vec', 1, equ_x_dim);
            equ_PI = zeros(1,5);
            for j=0:2
                for k=0:2
                    equ_PI(j+k+1) = equ_PI(j+k+1) + HMM_MODEL.PI(j+1)*HMM_MODEL.PI(k+1);
                end
            end
            %%            equ_PI = [HMM_MODEL.PI' HMM_MODEL.PI(1:end-1)']; equ_PI = equ_PI ./ sum(equ_PI);
            equ_SIGMA = HMM_MODEL.SIGMA*sqrt(2.0);   %% Not so good ... adding two noisy variables
            equ_N = ones(HMM_MODEL.x_dim*2-1, 1);
            equ_MEW = HMM_MODEL.MEW;  % Here we should sum the intensities of A and B - could be ambiguous
            equ_MEW(1) = 2*HMM_MODEL.MEW(1);
            equ_MEW(2) = HMM_MODEL.MEW(1)+HMM_MODEL.MEW(2);
            equ_MEW(3) = 0.5*( HMM_MODEL.MEW(1)+2*HMM_MODEL.MEW(2)+HMM_MODEL.MEW(3) );
            equ_MEW(4) = 0.5*( HMM_MODEL.MEW(1)+HMM_MODEL.MEW(2)+HMM_MODEL.MEW(3)+HMM_MODEL.MEW(4) );
            equ_MEW(5) = 0.5*( HMM_MODEL.MEW(1)+2*HMM_MODEL.MEW(3)+HMM_MODEL.MEW(5) );
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Learn model from data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(learn_flag)
            ttt_EM = cputime;
            [COMP_HMM_MODELS{i}.PI COMP_HMM_MODELS{i}.M COMP_HMM_MODELS{i}.N COMP_HMM_MODELS{i}.MEW COMP_HMM_MODELS{i}.SIGMA COMP_HMM_MODELS{i}.LogScore] = ...
                TrainHMMFromDataEMMatlab(OUT_VEC, OUT_VEC_B, LOC_VEC, use_locations, HMM_x_dim, HMM_y_dim, ...
                1, do_fold_change, mean_vec_rep, std_vec_rep+1, HMM_MODEL.PLACE_M(:,1:length(OUT_VEC)), ...
                HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.M_UPPERBOUNDS, HMM_MODEL.USE_BOUNDS, ...
                num_EM_iters, num_EM_starting_points, EM_tolerance);
            ttt_EM = cputime - ttt_EM
            COMP_HMM_MODELS{i}.M
            COMP_HMM_MODELS{i}.MEW
            COMP_HMM_MODELS{i}.SIGMA

            % Calculate Errors of the learned versus the correct model. Note : here we do not care WHERE is the error - very bad !
            PARAMS_ERROR.M_SQR_MEAN(:,:,points_ind) = PARAMS_ERROR.M_SQR_MEAN(:,:,points_ind) + ( (COMP_HMM_MODELS{i}.M-HMM_MODEL.M) ./ HMM_MODEL.M ) .^ 2;
            PARAMS_ERROR.MEW_SQR_MEAN(:,points_ind) = PARAMS_ERROR.MEW_SQR_MEAN(:,points_ind) + ( (COMP_HMM_MODELS{i}.MEW-HMM_MODEL.MEW)./MEW_GAP ) .^ 2;
            PARAMS_ERROR.SIGMA_SQR_MEAN(:,points_ind) = PARAMS_ERROR.SIGMA_SQR_MEAN(:,points_ind) + ( (COMP_HMM_MODELS{i}.SIGMA-HMM_MODEL.SIGMA') ./ HMM_MODEL.SIGMA' ) .^ 2;

            PARAMS_ERROR.M_SQR_STD(:,:,points_ind) = PARAMS_ERROR.M_SQR_STD(:,:,points_ind) + ( (COMP_HMM_MODELS{i}.M-HMM_MODEL.M) ./ HMM_MODEL.M ) .^ 4;
            PARAMS_ERROR.MEW_SQR_STD(:,points_ind) = PARAMS_ERROR.MEW_SQR_STD(:,points_ind) + ( COMP_HMM_MODELS{i}.MEW-HMM_MODEL.MEW ) .^ 4;
            PARAMS_ERROR.SIGMA_SQR_STD(:,points_ind) = PARAMS_ERROR.SIGMA_SQR_STD(:,points_ind) + ( (COMP_HMM_MODELS{i}.SIGMA-HMM_MODEL.SIGMA') ./ HMM_MODEL.SIGMA' ) .^ 4;

            % Calculate relative error distance
            KL_dist_seq_len=length(HMM_MODEL.PLACE_M)-1; % We must never allow more than length of PLACE_M here !!!!
            CUR_KL = ComputeHMMKLDistance(HMM_MODEL.SPECIAL_MODELS_FLAG,  KL_dist_seq_len, num_KL_dist_iters, ...
                max(HMM_MODEL.PLACE_M, HMM_MODEL.SPECIAL_MODELS_FLAG), ...
                HMM_MODEL.PI, HMM_MODEL.M, HMM_MODEL.N, ...
                HMM_MODEL.MEW, HMM_MODEL.SIGMA, ...
                COMP_HMM_MODELS{i}.PI, COMP_HMM_MODELS{i}.M, COMP_HMM_MODELS{i}.N, COMP_HMM_MODELS{i}.MEW, ...
                COMP_HMM_MODELS{i}.SIGMA)
            CLASS_ERROR.KL_MEAN(points_ind) =  CLASS_ERROR.KL_MEAN(points_ind) + CUR_KL;
            CLASS_ERROR.KL_STD(points_ind) =  CLASS_ERROR.KL_STD(points_ind) + (CUR_KL^2);

            % % %                         COMP_HMM_MODELS{i}.MEW=HMM_MODEL.MEW
        else
            COMP_HMM_MODELS{i} = CopyModels(HMM_MODEL);             % Here just copy the model parameters
        end  % if learn_flag

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Learn equivalent model from sum of copy number data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(learn_equ_flag)
            % Now also learn an 'equivalent' model, just to see that things are working fine
            [equ_MODEL.PI equ_MODEL.M equ_MODEL.N equ_MODEL.MEW equ_MODEL.SIGMA equ_MODEL.LogScore] = ...
                TrainHMMFromDataEMMatlab(OUT_VEC+OUT_VEC_B, OUT_VEC+OUT_VEC_B, LOC_VEC, use_locations, 2*HMM_x_dim-1, HMM_y_dim, ...
                0, do_fold_change, mean_vec_rep, std_vec_rep+1, HMM_MODEL.PLACE_M, ...
                0, HMM_MODEL.M_UPPERBOUNDS, HMM_MODEL.USE_BOUNDS, ...
                num_EM_iters, num_EM_starting_points, EM_tolerance);
            EQU_CUR_KL = ComputeHMMKLDistance(0,  KL_dist_seq_len, num_KL_dist_iters, ...
                HMM_MODEL.PLACE_M, ...
                equ_PI, equ_M, equ_N, ...
                equ_MEW, equ_SIGMA, ...
                equ_MODEL.PI, equ_MODEL.M, equ_MODEL.N, equ_MODEL.MEW, ...
                equ_MODEL.SIGMA)
            CLASS_ERROR.EQU_KL_MEAN(points_ind) =  CLASS_ERROR.EQU_KL_MEAN(points_ind) + EQU_CUR_KL;
            CLASS_ERROR.EQU_KL_STD(points_ind) =  CLASS_ERROR.EQU_KL_STD(points_ind) + (EQU_CUR_KL^2);

            % Calculate the 'distance' between the learned and original equivalent models
            PARAMS_ERROR.EQU_M_SQR_MEAN(:,:,points_ind) = PARAMS_ERROR.EQU_M_SQR_MEAN(:,:,points_ind) + ( (equ_MODEL.M-equ_M) ./ equ_M ) .^ 2;
            PARAMS_ERROR.EQU_MEW_SQR_MEAN(:,points_ind) = PARAMS_ERROR.EQU_MEW_SQR_MEAN(:,points_ind) + ( (equ_MODEL.MEW-equ_MEW)./MEW_GAP ) .^ 2;
            PARAMS_ERROR.EQU_SIGMA_SQR_MEAN(:,points_ind) = PARAMS_ERROR.EQU_SIGMA_SQR_MEAN(:,points_ind) + ( (equ_MODEL.SIGMA-equ_SIGMA') ./ equ_SIGMA' ) .^ 2;
            PARAMS_ERROR.EQU_M_SQR_STD(:,:,points_ind) = PARAMS_ERROR.EQU_M_SQR_STD(:,:,points_ind) + ( (equ_MODEL.M-equ_M) ./ equ_M ) .^ 4;
            PARAMS_ERROR.EQU_MEW_SQR_STD(:,points_ind) = PARAMS_ERROR.EQU_MEW_SQR_STD(:,points_ind) + ( equ_MODEL.MEW-equ_MEW ) .^ 4;
            PARAMS_ERROR.EQU_SIGMA_SQR_STD(:,points_ind) = PARAMS_ERROR.EQU_SIGMA_SQR_STD(:,points_ind) + ( (equ_MODEL.SIGMA-equ_SIGMA') ./ equ_SIGMA' ) .^ 4;

            equ_PI = equ_MODEL.PI;
            equ_M = equ_MODEL.M;
            equ_N = equ_MODEL.N;
            equ_MEW = equ_MODEL.MEW;
            equ_SIGMA = equ_MODEL.SIGMA;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now calculate Viterbi errors ('classification errors')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        mex FindBestPathViterbi.cpp hmm_chrom_funcs.cpp;
        [Viterbi_Path{i} Gamma_Probs{i}] = ...
            FindBestPathViterbi(OUT_VEC, OUT_VEC_B, LOC_VEC, ...
            use_locations, COMP_HMM_MODELS{i}.PI, COMP_HMM_MODELS{i}.M, COMP_HMM_MODELS{i}.N, ...
            COMP_HMM_MODELS{i}.MEW, COMP_HMM_MODELS{i}.SIGMA, max(HMM_MODEL.PLACE_FLAG, HMM_MODEL.SPECIAL_MODELS_FLAG), ...
            HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.SNP_specific, mean_vec_rep, std_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M);  % add the place variables all_chr_exp_arr(j,:)
%%         figure; imagesc(Gamma_Probs{2}(:,1:20)); colorbar;
        Viterbi_Path{i} = Viterbi_Path{i}'; %%Gamma_Probs{i} = Gamma_Probs{i}';
        % Get from the marginal Gamma probabilities, the argmax of the copy
        % numbers, geotypes etc.
        if( size(HMM_MODEL.PLACE_M,1) == 4 ) % Transpose if needed
            HMM_MODEL.PLACE_M = HMM_MODEL.PLACE_M';
        end
        if(HMM_MODEL.SPECIAL_MODELS_FLAG)
            [VitStruct ProbsStruct] = GetBestMarginalPredictions(Gamma_Probs{i}, HMM_MODEL, do_couples, joint_flag);
%%            [VitStruct_Singles ProbsStruct_Singels] = GetBestMarginalPredictions(Gamma_Probs{i}, HMM_MODEL, 0, joint_flag);            
            [Viterbi_Path_SINGLETONS{i} Gamma_Probs_SINGLETONS{i}] = ...
                FindBestPathViterbi(OUT_VEC, OUT_VEC_B, LOC_VEC, ...
                use_locations, COMP_HMM_MODELS{i}.PI, COMP_HMM_MODELS{i}.M, COMP_HMM_MODELS{i}.N, ...
                COMP_HMM_MODELS{i}.MEW, COMP_HMM_MODELS{i}.SIGMA, max(HMM_MODEL.PLACE_FLAG, HMM_MODEL.SPECIAL_MODELS_FLAG), ...
                HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.SNP_specific, mean_vec_rep, std_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M_SINGLETONS);  % add the place variables all_chr_exp_arr(j,:)
            Viterbi_Path_SINGLETONS{i} = Viterbi_Path_SINGLETONS{i}'; %%Gamma_Probs{i} = Gamma_Probs{i}';
            % Get from the marginal Gamma probabilities, the argmax of the copy
            % numbers, geotypes etc.
            %%       if( size(HMM_MODEL.PLACE_M,1) == 4 ) % Transpose if needed
            %%           HMM_MODEL.PLACE_M = HMM_MODEL.PLACE_M';
            %%       end
            [VitStruct_SINGLETONS ProbsStruct_SINGLETONS] = GetBestMarginalPredictions(Gamma_Probs_SINGLETONS{i}, HMM_MODEL, do_couples, joint_flag);

            [Viterbi_Path_IID{i} Gamma_Probs_IID{i}] = ...
                FindBestPathViterbi(OUT_VEC, OUT_VEC_B, LOC_VEC, ...
                use_locations, COMP_HMM_MODELS{i}.PI, COMP_HMM_MODELS{i}.M, COMP_HMM_MODELS{i}.N, ...
                COMP_HMM_MODELS{i}.MEW, COMP_HMM_MODELS{i}.SIGMA, max(HMM_MODEL.PLACE_FLAG, HMM_MODEL.SPECIAL_MODELS_FLAG), ...
                HMM_MODEL.SPECIAL_MODELS_FLAG, HMM_MODEL.SNP_specific, mean_vec_rep, std_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M_IID);  % add the place variables all_chr_exp_arr(j,:)
            Viterbi_Path_IID{i} = Viterbi_Path_IID{i}';
            % Get from the marginal Gamma probabilities, the argmax of the copy
            % numbers, geotypes etc.
            %%        if( size(HMM_MODEL.PLACE_M,1) == 4 ) % Transpose if needed
            %%            HMM_MODEL.PLACE_M = HMM_MODEL.PLACE_M';
            %%        end
            [VitStruct_IID ProbsStruct_IID] = GetBestMarginalPredictions(Gamma_Probs_IID{i}, HMM_MODEL, do_couples, joint_flag);

            % Here use the equivalent true model, rather than the learned one ..
            % Here we do not use genotype information and just look at total intensity
            [Viterbi_Path_equ{i} Gamma_Probs_equ{i}] = ...
                FindBestPathViterbi(OUT_VEC+OUT_VEC_B, OUT_VEC+OUT_VEC_B, LOC_VEC, ...
                use_locations, equ_PI, equ_M, equ_N, ...
                equ_MEW, equ_SIGMA, 0, ...
                0 , HMM_MODEL.SNP_specific, mean_vec_rep, std_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M);  % add the place variables all_chr_exp_arr(j,:);
            Viterbi_Path_equ{i} = Viterbi_Path_equ{i}'; %%Gamma_Probs{i} = Gamma_Probs{i}';

            V_VEC = VitStruct; %% GetMeaningfullCopyNumbersAndGenotypes(Viterbi_Path{i});% Transfer V_VEC to its 'real' values
            V_VEC.Gamma_Probs = Gamma_Probs{i}; % Take the gamma probs also ..
            V_VEC.copy_equ = Viterbi_Path_equ{i};
            V_VEC_SINGLETONS = VitStruct_SINGLETONS; %% GetMeaningfullCopyNumbersAndGenotypes(Viterbi_Path_SINGLETONS{i});% Transfer V_VEC to its 'real' values
            V_VEC_SINGLETONS.Gamma_Probs = Gamma_Probs_SINGLETONS{i}; % Take the gamma probs also ..
            V_VEC_SINGLETONS.copy_equ = Viterbi_Path_equ{i}; % WRONG! (taken from the full information model)
            V_VEC_IID = VitStruct_IID; %% GetMeaningfullCopyNumbersAndGenotypes(Viterbi_Path_IID{i});% Transfer V_VEC to its 'real' values
            V_VEC_IID.Gamma_Probs = Gamma_Probs_IID{i}; % Take the gamma probs also ..
            V_VEC_IID.copy_equ = Viterbi_Path_equ{i};  % WRONG! (taken from the full information model)

            CLASS_ERROR = ComputeClassificationError(X_VEC, V_VEC, HMM_MODEL, CLASS_ERROR, points_ind);        % Now find the classification error
            CLASS_ERROR_SINGLETONS = ComputeClassificationError(X_VEC, V_VEC_SINGLETONS, HMM_MODEL, CLASS_ERROR_SINGLETONS, points_ind);        % Now find the classification error
            CLASS_ERROR_IID = ComputeClassificationError(X_VEC, V_VEC_IID, HMM_MODEL, CLASS_ERROR_IID, points_ind);        % Now find the classification error

        end

        % Here the error is based on 0.5 threshold clipping !!!
        %       Gamma_Clips =  (Gamma_Probs{i}(2,:) > 0.5);
        %       CLASS_ERROR.BAYES_MEAN(points_ind) = CLASS_ERROR.BAYES_MEAN(points_ind) + (sum(abs(X_VEC'-Gamma_Clips)) / length(X_VEC));
        %       CLASS_ERROR.BAYES_STD(points_ind) = CLASS_ERROR.BAYES_STD(points_ind) + (sum(abs(X_VEC'-Gamma_Clips)) / length(X_VEC))^2;
        %         if( (i == num_iters) && (num_points == num_points_vec(end)) )
        %             figure; subplot(2,2,1); plot(X_VEC.alpha_copy, '.'); axis([0 500 0 2]); legend('\alpha copy');
        %             subplot(2,2,2); plot(X_VEC.beta_copy, '.'); axis([0 500 0 2]); legend('\beta copy');
        %             subplot(2,2,3); plot(V_VEC.alpha_copy, '.'); axis([0 500 0 2]); legend('recovered \alpha copy');
        %             subplot(2,2,4); plot(V_VEC.beta_copy, '.'); axis([0 500 0 2]); legend('recovered \beta copy');
        %
        %             figure; subplot(2,3,1); plot(min(X_VEC.alpha_copy, X_VEC.beta_copy), '.'); axis([0 500 0 2]); legend('min copy');
        %             subplot(2,3,2); plot(max(X_VEC.alpha_copy, X_VEC.beta_copy), '.'); axis([0 500 0 2]); legend('max copy');
        %             subplot(2,3,3); plot(X_VEC.alpha_copy+X_VEC.beta_copy, '.'); axis([0 500 0 4]); legend('sum copy');
        %             subplot(2,3,4); plot(min(V_VEC.alpha_copy, V_VEC.beta_copy), '.'); axis([0 500 0 2]); legend('recovered min copy');
        %             subplot(2,3,5); plot(max(V_VEC.alpha_copy, V_VEC.beta_copy), '.'); axis([0 500 0 2]); legend('recovered max copy');
        %             subplot(2,3,6); plot(V_VEC.alpha_copy+V_VEC.beta_copy, '.'); axis([0 500 0 4]); legend('recovered sum copy');
        %         end
    end % iter loop

    points_ind = points_ind + 1;
    done_num_points = num_points
end  % num_points loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing: Normalize the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS_ERROR = NormalizeErrors(CLASS_ERROR, num_iters);
CLASS_ERROR_SINGLETONS = NormalizeErrors(CLASS_ERROR_SINGLETONS, num_iters);
CLASS_ERROR_IID = NormalizeErrors(CLASS_ERROR_IID, num_iters);
PARAMS_ERROR = NormalizeParamsErrors(PARAMS_ERROR, num_iters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the plots showing performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorvec = 'bgrkmcxo:.';

% Plot classification error
if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
    figure; subplot(2,1,1); hold on; plot( (X_VEC.A_copy - V_VEC.A_copy).^2-0.05, '.'); title('Classification Error probability vectors');
    plot( (X_VEC.B_copy - V_VEC.B_copy).^2+0.05, 'r.'); plot( (V_VEC.total_copy-X_VEC.total_copy).^2+0.1,'g.');
    if(learn_equ_flag)
        plot( (V_VEC.copy_equ-X_VEC.copy).^2, 'm.');
        legend('A copy error', 'B copy error', 'Tot. copy error', 'Tot. copy error equ.');
    else
        legend('A copy error', 'B copy error', 'Tot. copy error');
    end
    xlabel('Time (t)'); ylabel('Square Error');
    subplot(2,1,2); hold on; plot( 0.02+(X_VEC.alpha_genotype - V_VEC.alpha_genotype).^2, '.');
    plot( (X_VEC.beta_genotype - V_VEC.beta_genotype).^2, 'r.');
    
    X_VEC_tot_geno = bitget(X_VEC.joint_genotype, 1) + bitget(X_VEC.joint_genotype, 2);
    V_VEC_tot_geno = bitget(V_VEC.joint_genotype, 1) + bitget(V_VEC.joint_genotype, 2);
    ERR_VEC_geno = abs(X_VEC_tot_geno - V_VEC_tot_geno);    
    ERR_VEC_geno = ERR_VEC_geno .* (1 - (X_VEC.total_copy == V_VEC.total_copy) .* (X_VEC.total_copy==1));
    ERR_VEC_geno = ERR_VEC_geno .* (1 - (X_VEC.total_copy == V_VEC.total_copy) .* (X_VEC.total_copy==0));
    % We cannot expect to get the correct genotype if one of the copy number is zero - so the convention is to put here -1
     
%%    zero_inds = intersect(find(X_VEC.A_copy == V_VEC.A_copy), find(X_VEC.B_copy == V_VEC.B_copy)); 
    zero_inds = intersect(find(X_VEC.alpha_copy == V_VEC.alpha_copy), find(X_VEC.beta_copy == V_VEC.beta_copy)); 
    zero_inds = intersect(zero_inds, find(min(X_VEC.A_copy, X_VEC.B_copy) == 0));
    ERR_VEC_geno(zero_inds) = -1;
    
    
    plot( 0.04 + ERR_VEC_geno, 'g.');
    legend('\alpha genotype error', '\beta genotype error', 'joint genotype error'); xlabel('Time (t)'); ylabel('Square Error');

    figure; hold on; title('Classification (Viterbi) Error Probability Vector A,B');
    xlabel('Num. Samples'); ylabel('CLASS. Error');
    errorbar(num_points_vec, CLASS_ERROR.A_MEAN, CLASS_ERROR.A_STD);
    errorbar(num_points_vec+1, CLASS_ERROR.B_MEAN, CLASS_ERROR.B_STD,'r');
    errorbar(num_points_vec+2, CLASS_ERROR.MEAN, CLASS_ERROR.STD,'g');
    if(learn_equ_flag)
        errorbar(num_points_vec+3, CLASS_ERROR.MEAN_equ, CLASS_ERROR.STD_equ,'m');
        legend('A copy', 'B copy', 'Total copy', 'Total copy equ.');
    else
        legend('A copy', 'B copy', 'Total copy'); 
    end
    
%     % Another Debugging Figures    
%     figure; hold on; plot(OUT_VEC, '.'); plot(OUT_VEC_B, 'r.'); 
%     plot(min(X_VEC.joint_genotype, 2)+0.03, 'g.');
%     legend('A Copy', 'B Copy', 'genotypes'); xlabel('Time (t)'); ylabel('Intensity');    
%     figure; hold on; title('A Vs. B Intensities'); xlabel('A Intensity'); ylabel('B Intensity');
%     plot(OUT_VEC(X_VEC.joint_genotype == 0), OUT_VEC_B(X_VEC.joint_genotype == 0), '.');
%     plot(OUT_VEC(X_VEC.joint_genotype == 1), OUT_VEC_B(X_VEC.joint_genotype == 1), 'r.');
%     plot(OUT_VEC(X_VEC.joint_genotype == 3), OUT_VEC_B(X_VEC.joint_genotype == 3), 'g.');
%     legend('AA', 'AB', 'BB');
%     fff = find(X_VEC.A_copy ~= V_VEC.A_copy);
%     figure; imagesc(Gamma_Probs{2}(:,fff)); colorbar; title('BAD 1<->2 Places');
%     figure; imagesc([A_copy_tab' B_copy_tab']); colorbar; title('A Copy              B Copy');    
%     % Now take just random ..
%     ggg = find( X_VEC.A_copy .*V_VEC.A_copy == 1);  
%     figure; imagesc(Gamma_Probs{2}(:,ggg)); colorbar; title('GOOD 1<->1 Places');
%     ggg = find( X_VEC.A_copy .*V_VEC.A_copy == 4);  
%     figure; imagesc(Gamma_Probs{2}(:,ggg)); colorbar; title('GOOD 2<->2 Places');    
%     [G_MAX G_MAX_IND] = max(Gamma_Probs{2});
%     unique(G_MAX_IND)
%     for iii=0:4
%         UNI{iii+1} = unique(G_MAX_IND(X_VEC.total_copy == iii))
%         UNI_MAT{iii+1} = find(Total_copy_tab == iii);
%     end
%     % Debug plotting - to be removed ..
%     [G_A_MAX G_A_MAX_IND] = max(ProbsStruct.A_copy');
%     [G_B_MAX G_B_MAX_IND] = max(ProbsStruct.B_copy');
%     [G_TOT_MAX G_TOT_MAX_IND] = max(ProbsStruct.total_copy');    
%     figure; plot(X_VEC.A_copy + 0.1*rand(length(X_VEC.A_copy),1), OUT_VEC + 0.1*rand(length(X_VEC.A_copy),1), '.'); title('Input dist');
%     figure; plot(X_VEC.A_copy + 0.1*rand(length(X_VEC.A_copy),1), A_copy_tab(G_MAX_IND)' + 0.1*rand(length(X_VEC.A_copy),1), '.'); title('Output dist');
%     figure; plot(X_VEC.A_copy + 0.1*rand(length(X_VEC.A_copy),1), G_TOT_MAX_IND' + 0.1*rand(length(X_VEC.A_copy),1), '.'); title('Output dist G TOT');
%     figure; plot(X_VEC.A_copy + 0.1*rand(length(X_VEC.A_copy),1), G_A_MAX_IND' + 0.1*rand(length(X_VEC.A_copy),1), '.'); title('Output dist G A');
%     figure; plot(X_VEC.B_copy + 0.1*rand(length(X_VEC.B_copy),1), G_B_MAX_IND' + 0.1*rand(length(X_VEC.B_copy),1), '.'); title('Output dist G B');        
%     figure; plot(X_VEC.A_copy + 0.1*rand(length(X_VEC.A_copy),1), V_VEC.A_copy + 0.1*rand(length(X_VEC.A_copy),1), '.'); title('Output dist Again');
%     figure; plot(X_VEC.A_copy + 0.1*rand(length(X_VEC.A_copy),1), VitStruct.A_copy + 0.1*rand(length(X_VEC.A_copy),1), '.'); title('Output dist Again Singles');

    PlotClassErrorComparison(CLASS_ERROR, CLASS_ERROR_SINGLETONS, CLASS_ERROR_IID, num_points_vec);
else
    % Plot the Classification (Viterbi) error
    figure; subplot(2,1,1); hold on; title('Classification (Viterbi) Error Probability Vector'); xlabel('Num. Samples'); ylabel('CLASS. Error');
    errorbar(num_points_vec, CLASS_ERROR.MEAN, CLASS_ERROR.STD);
    % Plot the Bayesian Classification (forward) error
    subplot(2,1,2); hold on; title('Classification (forward) "Error Probability" Vector'); xlabel('Num. Samples'); ylabel('forward CLASS. Error');
    errorbar(num_points_vec, CLASS_ERROR.BAYES_MEAN, CLASS_ERROR.BAYES_STD);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot learned model parameters for showing performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(learn_flag)
    % Plot M errors
    figure; subplot(2,2,1); hold on;  title('Relative Errors in M - Transition Matrix'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    for i=1:HMM_x_dim
        for j=1:HMM_x_dim
            errorbar(num_points_vec, reshape(PARAMS_ERROR.M_SQR_MEAN(i,j,:), 1, K), reshape(PARAMS_ERROR.M_SQR_STD(i,j,:), 1, K), colorvec((i-1)*HMM_x_dim+j));
        end
    end
    subplot(2,2,2); hold on; title('Relative Errors in MEW (divided by mew gaps) - Expectation Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    % Plot MEW errors
    for i=1:HMM_x_dim
        errorbar(num_points_vec, PARAMS_ERROR.MEW_SQR_MEAN(i,:), PARAMS_ERROR.MEW_SQR_STD(i,:), colorvec(i));
    end
    subplot(2,2,3); hold on; title('Relative Errors in SIGMA - Standard Error Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    % Plot MEW errors
    for i=1:HMM_x_dim
        errorbar(num_points_vec, PARAMS_ERROR.SIGMA_SQR_MEAN(i,:), PARAMS_ERROR.SIGMA_SQR_STD(i,:), colorvec(i));
    end
    % Plot the KL distance
    subplot(2,2,4); hold on; title('Kullback-Leibler Error Vector'); xlabel('Num. Samples'); ylabel('KL Error');
    errorbar(num_points_vec, CLASS_ERROR.KL_MEAN, CLASS_ERROR.KL_STD);
end

if(learn_equ_flag)
    % Plot the equivalent models:
    % Plot equ M errors
    figure; subplot(2,2,1); hold on;  title('EQU Relative Errors in M - Transition Matrix'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    for i=1:HMM_x_dim
        for j=1:HMM_x_dim
            errorbar(num_points_vec, reshape(PARAMS_ERROR.EQU_M_SQR_MEAN(i,j,:), 1, K), reshape(PARAMS_ERROR.EQU_M_SQR_STD(i,j,:), 1, K), colorvec((i-1)*HMM_x_dim+j));
        end
    end
    subplot(2,2,2); hold on; title('EQU Relative Errors in MEW (divided by mew gaps) - Expectation Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    % Plot equ MEW errors
    for i=1:HMM_x_dim
        errorbar(num_points_vec, PARAMS_ERROR.EQU_MEW_SQR_MEAN(i,:), PARAMS_ERROR.EQU_MEW_SQR_STD(i,:), colorvec(i));
    end
    subplot(2,2,3); hold on; title('EQU Relative Errors in SIGMA - Standard Error Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    % Plot equ SIGMA errors
    for i=1:HMM_x_dim
        errorbar(num_points_vec, PARAMS_ERROR.EQU_SIGMA_SQR_MEAN(i,:), PARAMS_ERROR.EQU_SIGMA_SQR_STD(i,:), colorvec(i));
    end
    % Plot equ KL distance
    subplot(2,2,4); hold on; title('EQU Kullback-Leibler Error Vector'); xlabel('Num. Samples'); ylabel('KL Error');
    errorbar(num_points_vec, CLASS_ERROR.EQU_KL_MEAN, CLASS_ERROR.EQU_KL_STD);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxillary functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the correct copy number and genotypes from a 'packed' joint form
function  X_VEC = GetMeaningfullCopyNumbersAndGenotypes(all_data);
% Convension: A - 0, B - 1

X_VEC.all_data = all_data;
X_VEC.alpha_genotype = bitget(all_data,1); X_VEC.beta_genotype = bitget(all_data,5);
X_VEC.joint_genotype = X_VEC.alpha_genotype+X_VEC.beta_genotype + X_VEC.alpha_genotype.*X_VEC.beta_genotype;
X_VEC.alpha_copy = bitget(all_data,2) + 2*bitget(all_data,3) + 4*bitget(all_data,4);
X_VEC.beta_copy = bitget(all_data,6) + 2*bitget(all_data,7) + 4*bitget(all_data,8);
X_VEC.total_copy = X_VEC.alpha_copy+X_VEC.beta_copy;
X_VEC.A_copy = (1-X_VEC.alpha_genotype).*X_VEC.alpha_copy + (1-X_VEC.beta_genotype).*X_VEC.beta_copy;
X_VEC.B_copy = X_VEC.alpha_genotype.*X_VEC.alpha_copy + X_VEC.beta_genotype.*X_VEC.beta_copy;



% Update the values of the classifcation error based on the current values of the two vectors
function CLASS_ERROR = ComputeClassificationError(X_VEC, V_VEC, HMM_MODEL, CLASS_ERROR, points_ind)
% Now find the classification error
if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
    CLASS_ERROR.MEAN(points_ind) = CLASS_ERROR.MEAN(points_ind) + sum( (V_VEC.total_copy - X_VEC.total_copy).^2 ./ length(X_VEC.total_copy) );
    CLASS_ERROR.STD(points_ind) = CLASS_ERROR.STD(points_ind) + (sum( (V_VEC.total_copy - X_VEC.total_copy).^2 / length(X_VEC.total_copy) ))^2;
    CLASS_ERROR.MEAN_equ(points_ind) = CLASS_ERROR.MEAN_equ(points_ind) + sum( (V_VEC.copy_equ - X_VEC.total_copy).^2 ./ length(X_VEC.total_copy) );
    CLASS_ERROR.STD_equ(points_ind) = CLASS_ERROR.STD_equ(points_ind) + (sum( (V_VEC.copy_equ - X_VEC.total_copy).^2 / length(X_VEC.total_copy) ))^2;

    CLASS_ERROR.GENOTYPE_MEAN(points_ind) = CLASS_ERROR.GENOTYPE_MEAN(points_ind) + ...
        sum( ( (V_VEC.joint_genotype ~= X_VEC.joint_genotype)  .* (X_VEC.total_copy + V_VEC.total_copy > 0) ) ./ length(X_VEC.joint_genotype) );
    CLASS_ERROR.GENOTYPE_STD(points_ind) = CLASS_ERROR.GENOTYPE_STD(points_ind) + ...
        sum( ( (V_VEC.joint_genotype ~= X_VEC.joint_genotype) .* (X_VEC.total_copy + V_VEC.total_copy > 0)) ./ length(X_VEC.joint_genotype) )^2;

else
    CLASS_ERROR.MEAN(points_ind) = CLASS_ERROR.MEAN(points_ind) + (sum(Viterbi_Path{i} ~= X_VEC) / length(X_VEC));
    CLASS_ERROR.STD(points_ind) = CLASS_ERROR.STD(points_ind) + (sum(Viterbi_Path{i} ~= X_VEC) / length(X_VEC))^2;
end
CLASS_ERROR.A_MEAN(points_ind) = CLASS_ERROR.A_MEAN(points_ind) + sum( (X_VEC.A_copy-V_VEC.A_copy).^2 ./ length(X_VEC.A_copy) );
CLASS_ERROR.B_MEAN(points_ind) = CLASS_ERROR.B_MEAN(points_ind) + sum( (X_VEC.B_copy-V_VEC.B_copy).^2 ./ length(X_VEC.B_copy) );
CLASS_ERROR.A_STD(points_ind) = CLASS_ERROR.A_STD(points_ind) + (sum( (X_VEC.A_copy-V_VEC.A_copy).^2 ./ length(X_VEC.A_copy) ))^2;
CLASS_ERROR.B_STD(points_ind) = CLASS_ERROR.B_STD(points_ind) + (sum( (X_VEC.B_copy-V_VEC.B_copy).^2 ./ length(X_VEC.B_copy) ))^2;
% Here the error is weighted according to the probability given
CLASS_ERROR.BAYES_MEAN(points_ind) = CLASS_ERROR.BAYES_MEAN(points_ind) + (sum(abs(X_VEC.all_data'-V_VEC.Gamma_Probs(2,:))) / length(X_VEC.all_data));
CLASS_ERROR.BAYES_STD(points_ind) = CLASS_ERROR.BAYES_STD(points_ind) + (sum(abs(X_VEC.all_data'-V_VEC.Gamma_Probs(2,:))) / length(X_VEC.all_data))^2;


% Normalize the errors - just divide by the number of iterations
function CLASS_ERROR = NormalizeErrors(CLASS_ERROR, num_iters)
CLASS_ERROR.MEAN = CLASS_ERROR.MEAN ./ num_iters;
CLASS_ERROR.BAYES_MEAN = CLASS_ERROR.BAYES_MEAN ./ num_iters;
CLASS_ERROR.A_MEAN = CLASS_ERROR.A_MEAN ./ num_iters;
CLASS_ERROR.B_MEAN = CLASS_ERROR.B_MEAN ./ num_iters;
CLASS_ERROR.MEAN_equ = CLASS_ERROR.MEAN_equ ./ num_iters;
CLASS_ERROR.GENOTYPE_MEAN = CLASS_ERROR.GENOTYPE_MEAN ./ num_iters;

CLASS_ERROR.STD = sqrt(CLASS_ERROR.STD ./ num_iters - (CLASS_ERROR.MEAN .^ 2));
CLASS_ERROR.BAYES_STD = sqrt(CLASS_ERROR.BAYES_STD ./ num_iters - (CLASS_ERROR.BAYES_MEAN .^ 2));
CLASS_ERROR.A_STD = sqrt(CLASS_ERROR.A_STD ./ num_iters - (CLASS_ERROR.A_MEAN .^ 2));
CLASS_ERROR.B_STD = sqrt(CLASS_ERROR.B_STD ./ num_iters - (CLASS_ERROR.B_MEAN .^ 2));
CLASS_ERROR.STD_equ = sqrt(CLASS_ERROR.STD_equ ./ num_iters - (CLASS_ERROR.MEAN_equ .^ 2));
CLASS_ERROR.GENOTYPE_STD = sqrt(CLASS_ERROR.GENOTYPE_STD ./ num_iters - (CLASS_ERROR.GENOTYPE_MEAN .^ 2));

CLASS_ERROR.KL_MEAN = CLASS_ERROR.KL_MEAN ./ num_iters;
CLASS_ERROR.EQU_KL_MEAN = CLASS_ERROR.EQU_KL_MEAN ./ num_iters;
CLASS_ERROR.KL_STD = sqrt(CLASS_ERROR.KL_STD ./ num_iters - (CLASS_ERROR.KL_MEAN .^ 2));
CLASS_ERROR.EQU_KL_STD = sqrt(CLASS_ERROR.EQU_KL_STD ./ num_iters - (CLASS_ERROR.EQU_KL_MEAN .^ 2));

% Normalize the errors for parameters - just divide by the number of iterations
function PARAMS_ERROR = NormalizeParamsErrors(PARAMS_ERROR, num_iters)
PARAMS_ERROR.M_SQR_MEAN = PARAMS_ERROR.M_SQR_MEAN ./ num_iters;
PARAMS_ERROR.MEW_SQR_MEAN = PARAMS_ERROR.MEW_SQR_MEAN ./ num_iters;
PARAMS_ERROR.SIGMA_SQR_MEAN = PARAMS_ERROR.SIGMA_SQR_MEAN ./ num_iters;

PARAMS_ERROR.EQU_M_SQR_MEAN = PARAMS_ERROR.EQU_M_SQR_MEAN ./ num_iters;
PARAMS_ERROR.EQU_MEW_SQR_MEAN = PARAMS_ERROR.EQU_MEW_SQR_MEAN ./ num_iters;
PARAMS_ERROR.EQU_SIGMA_SQR_MEAN = PARAMS_ERROR.EQU_SIGMA_SQR_MEAN ./ num_iters;

PARAMS_ERROR.M_SQR_STD = sqrt(PARAMS_ERROR.M_SQR_STD ./ num_iters - (PARAMS_ERROR.M_SQR_MEAN .^ 2));
PARAMS_ERROR.MEW_SQR_STD = sqrt(PARAMS_ERROR.MEW_SQR_STD ./ num_iters - (PARAMS_ERROR.MEW_SQR_MEAN .^ 2));
PARAMS_ERROR.SIGMA_SQR_STD = sqrt(PARAMS_ERROR.SIGMA_SQR_STD ./ num_iters - (PARAMS_ERROR.SIGMA_SQR_MEAN .^ 2));

PARAMS_ERROR.EQU_M_SQR_STD = sqrt(PARAMS_ERROR.EQU_M_SQR_STD ./ num_iters - (PARAMS_ERROR.EQU_M_SQR_MEAN .^ 2));
PARAMS_ERROR.EQU_MEW_SQR_STD = sqrt(PARAMS_ERROR.EQU_MEW_SQR_STD ./ num_iters - (PARAMS_ERROR.EQU_MEW_SQR_MEAN .^ 2));
PARAMS_ERROR.EQU_SIGMA_SQR_STD = sqrt(PARAMS_ERROR.EQU_SIGMA_SQR_STD ./ num_iters - (PARAMS_ERROR.EQU_SIGMA_SQR_MEAN .^ 2));








% Initilize classification errors: Set everything to zero
function CLASS_ERROR = InitClassError(K)

CLASS_ERROR.MEAN = zeros(1, K);
CLASS_ERROR.BAYES_MEAN = zeros(1, K);
CLASS_ERROR.A_MEAN = zeros(1, K);
CLASS_ERROR.B_MEAN = zeros(1, K);
CLASS_ERROR.MEAN_equ = zeros(1, K);
CLASS_ERROR.GENOTYPE_MEAN = zeros(1, K);

CLASS_ERROR.STD = zeros(1, K);
CLASS_ERROR.BAYES_STD = zeros(1, K);
CLASS_ERROR.A_STD = zeros(1, K);
CLASS_ERROR.B_STD = zeros(1, K);
CLASS_ERROR.STD_equ = zeros(1, K);
CLASS_ERROR.GENOTYPE_STD = zeros(1, K);

CLASS_ERROR.KL_MEAN = zeros(1, K);
CLASS_ERROR.EQU_KL_MEAN = zeros(1, K);
CLASS_ERROR.KL_STD = zeros(1, K);
CLASS_ERROR.EQU_KL_STD = zeros(1, K);


% Initilize parameters errors: Set everything to zero
function PARAMS_ERROR = InitParamsError(K, HMM_MODEL)
PARAMS_ERROR.M_SQR_MEAN = zeros(HMM_MODEL.x_dim, HMM_MODEL.x_dim, K);
if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
    PARAMS_ERROR.MEW_SQR_MEAN = zeros(2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.SIGMA_SQR_MEAN = zeros(2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.EQU_M_SQR_MEAN = zeros(2*HMM_MODEL.x_dim-1, 2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.EQU_MEW_SQR_MEAN = zeros(2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.EQU_SIGMA_SQR_MEAN = zeros(2*HMM_MODEL.x_dim-1, K);
else
    PARAMS_ERROR.MEW_SQR_MEAN = zeros(HMM_MODEL.x_dim, K);
    PARAMS_ERROR.SIGMA_SQR_MEAN = zeros(HMM_MODEL.x_dim, K);
end

PARAMS_ERROR.M_SQR_STD = zeros(HMM_MODEL.x_dim, HMM_MODEL.x_dim, K);
if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
    PARAMS_ERROR.MEW_SQR_STD = zeros(2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.SIGMA_SQR_STD = zeros(2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.EQU_M_SQR_STD = zeros(2*HMM_MODEL.x_dim-1, 2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.EQU_MEW_SQR_STD = zeros(2*HMM_MODEL.x_dim-1, K);
    PARAMS_ERROR.EQU_SIGMA_SQR_STD = zeros(2*HMM_MODEL.x_dim-1, K);
else
    PARAMS_ERROR.MEW_SQR_STD = zeros(HMM_MODEL.x_dim, K);
    PARAMS_ERROR.SIGMA_SQR_STD = zeros(HMM_MODEL.x_dim, K);
end


% Here just copy the model parameters
function MODEL_TAR = CopyModels(MODEL_SRC);
MODEL_TAR = [];
MODEL_TAR.PI = MODEL_SRC.PI;
MODEL_TAR.M = MODEL_SRC.M;
MODEL_TAR.N = MODEL_SRC.N;
MODEL_TAR.PLACE_M = MODEL_SRC.PLACE_M;
MODEL_TAR.SPECIAL_MODELS_FLAG = MODEL_SRC.SPECIAL_MODELS_FLAG;
MODEL_TAR.MEW = MODEL_SRC.MEW;
MODEL_TAR.SIGMA = MODEL_SRC.SIGMA;


% Plot comparison: both copy number and genotypes
function  PlotClassErrorComparison(CLASS_ERROR, CLASS_ERROR_SINGLETONS, CLASS_ERROR_IID, num_points_vec);

% First plot the copy number errors
figure; subplot(2,2,1); hold on; title('Classification (Viterbi) Error Probability Vector A Copy');
xlabel('Num. Samples'); ylabel('CLASS. Error');
errorbar(num_points_vec, CLASS_ERROR.A_MEAN, CLASS_ERROR.A_STD);
errorbar(num_points_vec+1, CLASS_ERROR_SINGLETONS.A_MEAN, CLASS_ERROR_SINGLETONS.A_STD,'r');
errorbar(num_points_vec+2, CLASS_ERROR_IID.A_MEAN, CLASS_ERROR_IID.A_STD,'g');
legend('All Info.', 'Singletons', 'I.I.D.');

subplot(2,2,2); hold on; title('Classification (Viterbi) Error Probability Vector B Copy');
xlabel('Num. Samples'); ylabel('CLASS. Error');
errorbar(num_points_vec, CLASS_ERROR.B_MEAN, CLASS_ERROR.B_STD);
errorbar(num_points_vec+1, CLASS_ERROR_SINGLETONS.B_MEAN, CLASS_ERROR_SINGLETONS.B_STD,'r');
errorbar(num_points_vec+2, CLASS_ERROR_IID.B_MEAN, CLASS_ERROR_IID.B_STD,'g');
legend('All Info.', 'Singletons', 'I.I.D.');

subplot(2,2,3); hold on; title('Classification (Viterbi) Error Probability Vector Total Copy');
xlabel('Num. Samples'); ylabel('CLASS. Error');
errorbar(num_points_vec, CLASS_ERROR.MEAN, CLASS_ERROR.STD);
errorbar(num_points_vec+1, CLASS_ERROR_SINGLETONS.MEAN, CLASS_ERROR_SINGLETONS.STD,'r');
errorbar(num_points_vec+2, CLASS_ERROR_IID.MEAN, CLASS_ERROR_IID.STD,'g');
legend('All Info.', 'Singletons', 'I.I.D.');

% Now plot the genotypes error
subplot(2,2,4); hold on; title('Classification (Viterbi) Joint Genotype Error Probability Vector');
xlabel('Num. Samples'); ylabel('CLASS. Error');
errorbar(num_points_vec, CLASS_ERROR.GENOTYPE_MEAN, CLASS_ERROR.GENOTYPE_STD);
errorbar(num_points_vec+1, CLASS_ERROR_SINGLETONS.GENOTYPE_MEAN, CLASS_ERROR_SINGLETONS.GENOTYPE_STD,'r');
errorbar(num_points_vec+2, CLASS_ERROR_IID.GENOTYPE_MEAN, CLASS_ERROR_IID.GENOTYPE_STD,'g');
legend('All Info.', 'Singletons', 'I.I.D.');

