% Find what happens if we generate data from a given model and try to learn it again.
function [COMP_HMM_MODELS] = DetermineModelStability( HMM_MODEL, num_points_vec, num_iters, learn_flag)



% Make an 'equivalent' model which holds information only on the copy
% number
if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
    eq_x_dim = HMM_MODEL.x_dim*2-1; % Allow to add the A and B copy numbers
    eq_M = zeros(eq_x_dim);  % Build the 'equivalence' M markov matrix 
    prod_M = zeros(HMM_MODEL.x_dim.^2); 
    for i1=1:HMM_MODEL.x_dim
        for i2=1:HMM_MODEL.x_dim       
            for j1=1:HMM_MODEL.x_dim
                for j2=1:HMM_MODEL.x_dim
                    prod_M(i1+HMM_MODEL.x_dim*i2,j1+HMM_MODEL.x_dim*j2) = HMM_MODEL.M(i1,j1)* HMM_MODEL.M(i2,j2);
                    eq_M(i1+i2-1,j1+j2-1) = eq_M(i1+i2-1,j1+j2-1) + prod_M(i1+HMM_MODEL.x_dim*i2,j1+HMM_MODEL.x_dim*j2);
                end
            end
        end
    end
    eq_PI = [HMM_MODEL.PI' HMM_MODEL.PI(1:end-1)']; eq_PI = eq_PI ./ sum(eq_PI.^2);  
    eq_SIGMA = [HMM_MODEL.SIGMA HMM_MODEL.SIGMA(1:end-1)];
    eq_N = ones(HMM_MODEL.x_dim*2-1, 1); 
    eq_MEW = [0:HMM_MODEL.x_dim*2-2];
end


% Do many simulations, and each time determine the best learned model
COMP_HMM_MODELS = {};

% These structures contain the errors
M_SQR_ERROR_MEAN = zeros(length(HMM_MODEL.M), length(HMM_MODEL.M), length(num_points_vec));
MEW_SQR_ERROR_MEAN = zeros(length(HMM_MODEL.M), length(num_points_vec));
SIGMA_SQR_ERROR_MEAN = zeros(length(HMM_MODEL.M), length(num_points_vec));
KL_ERROR_MEAN = zeros(1, length(num_points_vec));
CLASS_ERROR_MEAN = zeros(1, length(num_points_vec));
BAYES_CLASS_ERROR_MEAN = zeros(1, length(num_points_vec));
CLASS_ERROR_A_MEAN = zeros(1, length(num_points_vec));
CLASS_ERROR_B_MEAN = zeros(1, length(num_points_vec));
CLASS_ERROR_MEAN_eq = zeros(1, length(num_points_vec));


M_SQR_ERROR_STD = zeros(length(HMM_MODEL.M), length(HMM_MODEL.M), length(num_points_vec));
MEW_SQR_ERROR_STD = zeros(length(HMM_MODEL.M), length(num_points_vec));
SIGMA_SQR_ERROR_STD = zeros(length(HMM_MODEL.M), length(num_points_vec));
KL_ERROR_STD = zeros(1, length(num_points_vec));
CLASS_ERROR_STD = zeros(1, length(num_points_vec));
BAYES_CLASS_ERROR_STD = zeros(1, length(num_points_vec));
CLASS_ERROR_A_STD = zeros(1, length(num_points_vec));
CLASS_ERROR_B_STD = zeros(1, length(num_points_vec));
CLASS_ERROR_STD_eq = zeros(1, length(num_points_vec));

% Parameters for KL distance
KL_dist_seq_len = 50000;
num_KL_dist_iters = 10;

points_ind = 1;
MEW_GAP  = max(HMM_MODEL.MEW) - min(HMM_MODEL.MEW); % The gap between the mew's for 'normalizing' the results

for num_points = num_points_vec
    for i=1:num_iters
        % Simulate sequence
        SAVE_M = HMM_MODEL.M; HMM_MODEL.M = eye(3); % Don't allow ANY transitions !!! 
        HMM_MODEL.PI = [0 1 0]; % Force the level to be one!!! 
        [X_VEC OUT_VEC OUT_VEC2] = SimulateSequenceFromModel(HMM_MODEL.PI, HMM_MODEL.M, HMM_MODEL.N, ...
            HMM_MODEL.MEW, HMM_MODEL.SIGMA, HMM_MODEL.PLACE_FLAG, ...
            HMM_MODEL.PLACE_M, HMM_MODEL.SPECIAL_MODELS_FLAG, num_points);
        % Change them back for the learning part
        HMM_MODEL.y_dim = 1;
        HMM_MODEL.Y_TYPE = 1;
        HMM_MODEL.M = SAVE_M;
        
        % Transfer X_VEC to its 'real' values
        X_VEC_alpha_genotype = bitget(X_VEC,1);
        X_VEC_beta_genotype = bitget(X_VEC,5);
        X_VEC_alpha_copynumber = bitget(X_VEC,2) + 2*bitget(X_VEC,3) + 4*bitget(X_VEC,4);
        X_VEC_beta_copynumber = bitget(X_VEC,6) + 2*bitget(X_VEC,7) + 4*bitget(X_VEC,8);
        X_VEC_copy = X_VEC_alpha_copynumber+X_VEC_beta_copynumber;
        X_VEC_A = X_VEC_alpha_genotype.*X_VEC_alpha_copynumber + ...
            X_VEC_beta_genotype.*X_VEC_beta_copynumber;
        X_VEC_B = (1-X_VEC_alpha_genotype).*X_VEC_alpha_copynumber + ...
            (1-X_VEC_beta_genotype).*X_VEC_beta_copynumber;


        % % %         % Now plot the simulation results
        % % %         figure; subplot(2,1,1); hold on;
        % % %         plot(X_VEC_alpha_genotype.*X_VEC_alpha_copynumber+0.06, 'g+');
        % % %         plot((1-X_VEC_alpha_genotype).*X_VEC_alpha_copynumber+0.06, 'go');
        % % %         plot(X_VEC_beta_genotype.*X_VEC_beta_copynumber-0.06, 'y+');
        % % %         plot((1-X_VEC_beta_genotype).*X_VEC_beta_copynumber-0.06, 'yo');
        % % %
        % % %         plot(X_VEC_alpha_copynumber-0.03, 'b');  plot(X_VEC_beta_copynumber+0.03, 'r');
        % % %         plot(X_VEC_alpha_copynumber+X_VEC_beta_copynumber+0.05, 'm');
        % % %  %%       plot(X_VEC_A+0.2, 'ko'); plot(X_VEC_B+0.2, 'cx');
        % % % %%%        plot(OUT_VEC, 'g.');  plot(OUT_VEC2, 'y.');
        % % %
        % % %         legend('alpha copy number', 'beta copy number', 'SUM', 'A genotype', 'B genotpye', 'A Intensity', 'B Intensity');
        % % %         xlabel('SNP location'); ylabel('Values');



        % Set all parameters
        LOC_VEC = 1:length(OUT_VEC);
        use_locations = 0;
        HMM_y_dim = 1;
        HMM_x_dim = length(HMM_MODEL.PI);
        do_fold_change = 0; % Note : Fold Change was all ready done !!!!
        mean_vec_rep = zeros(length(OUT_VEC),HMM_x_dim);
        std_vec_rep = zeros(length(OUT_VEC),HMM_x_dim);
        num_EM_iters = 100;
        num_EM_starting_points = 10;
        EM_tolerance = 0.000000001;

        % Learn model
        if(learn_flag)
            use_bounds = HMM_MODEL.USE_BOUNDS;
            [COMP_HMM_MODELS{i}.PI COMP_HMM_MODELS{i}.M COMP_HMM_MODELS{i}.N COMP_HMM_MODELS{i}.MEW COMP_HMM_MODELS{i}.SIGMA COMP_HMM_MODELS{i}.LogScore] = ...
                TrainHMMFromDataEMMatlab(OUT_VEC, LOC_VEC, use_locations, HMM_x_dim, HMM_y_dim, ...
                0, do_fold_change, mean_vec_rep, std_vec_rep, ...
                HMM_MODEL.M_UPPERBOUNDS, HMM_MODEL.USE_BOUNDS, ...
                num_EM_iters, num_EM_starting_points, EM_tolerance);

            % Calculate Errors of the learned versus the correct model. Note : here we do not care WHERE is the error - very bad !
            M_SQR_ERROR_MEAN(:,:,points_ind) = M_SQR_ERROR_MEAN(:,:,points_ind) + ( (COMP_HMM_MODELS{i}.M-HMM_MODEL.M) ./ HMM_MODEL.M ) .^ 2;
            MEW_SQR_ERROR_MEAN(:,points_ind) = MEW_SQR_ERROR_MEAN(:,points_ind) + ( (COMP_HMM_MODELS{i}.MEW-HMM_MODEL.MEW)/MEW_GAP ) .^ 2;
            SIGMA_SQR_ERROR_MEAN(:,points_ind) = SIGMA_SQR_ERROR_MEAN(:,points_ind) + ( (COMP_HMM_MODELS{i}.SIGMA-HMM_MODEL.SIGMA) ./ HMM_MODEL.SIGMA ) .^ 2;

            M_SQR_ERROR_STD(:,:,points_ind) = M_SQR_ERROR_STD(:,:,points_ind) + ( (COMP_HMM_MODELS{i}.M-HMM_MODEL.M) ./ HMM_MODEL.M ) .^ 4;
            MEW_SQR_ERROR_STD(:,points_ind) = MEW_SQR_ERROR_STD(:,points_ind) + ( COMP_HMM_MODELS{i}.MEW-HMM_MODEL.MEW ) .^ 4;
            SIGMA_SQR_ERROR_STD(:,points_ind) = SIGMA_SQR_ERROR_STD(:,points_ind) + ( (COMP_HMM_MODELS{i}.SIGMA-HMM_MODEL.SIGMA) ./ HMM_MODEL.SIGMA ) .^ 4;


            % Calculate relative error distance
            CUR_KL = ComputeHMMKLDistance(HMM_MODEL.PI, HMM_MODEL.M, HMM_MODEL.N, ...
                HMM_MODEL.MEW, HMM_MODEL.SIGMA, ...
                COMP_HMM_MODELS{i}.PI, COMP_HMM_MODELS{i}.M, COMP_HMM_MODELS{i}.N, COMP_HMM_MODELS{i}.MEW, ...
                COMP_HMM_MODELS{i}.SIGMA, KL_dist_seq_len, num_KL_dist_iters);


            KL_ERROR_MEAN(points_ind) =  KL_ERROR_MEAN(points_ind) + CUR_KL;
            KL_ERROR_STD(points_ind) =  KL_ERROR_STD(points_ind) + (CUR_KL^2);



            % Now calculate Viterbi errors ('classification errors')
            [Viterbi_Path{i} Gamma_Probs{i}] = ...
                FindBestPathViterbi(OUT_VEC, OUT_VEC2, LOC_VEC, ...
                use_locations, COMP_HMM_MODELS{i}.PI, COMP_HMM_MODELS{i}.M, COMP_HMM_MODELS{i}.N, ...
                COMP_HMM_MODELS{i}.MEW, COMP_HMM_MODELS{i}.SIGMA, HMM_MODEL.PLACE_FLAG, ...
                HMM_MODEL.SPECIAL_MODELS_FLAG, mean_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M);  % add the place variables all_chr_exp_arr(j,:)
            Viterbi_Path{i} = Viterbi_Path{i}'; %%Gamma_Probs{i} = Gamma_Probs{i}';
        else

            %%      doing_Viterbi = 999

            % Here use the real true model, rather than the learned one ..
            [Viterbi_Path{i} Gamma_Probs{i}] = ...
                [VVVV GGGG] = ....
                FindBestPathViterbi(OUT_VEC, OUT_VEC2, LOC_VEC, ...
                use_locations, HMM_MODEL.PI, HMM_MODEL.M, HMM_MODEL.N, ...
                HMM_MODEL.MEW, HMM_MODEL.SIGMA, HMM_MODEL.PLACE_FLAG, ...
                HMM_MODEL.SPECIAL_MODELS_FLAG , mean_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M);  % add the place variables all_chr_exp_arr(j,:);
            Viterbi_Path{i} = Viterbi_Path{i}'; %%Gamma_Probs{i} = Gamma_Probs{i}';

        
            % Here use the equivalent true model, rather than the learned one ..
            % Here do not use genotype information and just look at total intensity
            [Viterbi_Path_eq{i} Gamma_Probs_eq{i}] = ...
                FindBestPathViterbi(OUT_VEC+OUT_VEC2, OUT_VEC+OUT_VEC2, LOC_VEC, ...
                use_locations, eq_PI, eq_M, eq_N, ...
                eq_MEW, eq_SIGMA, 0, ...
                0 , mean_vec_rep, std_vec_rep, HMM_MODEL.PLACE_M);  % add the place variables all_chr_exp_arr(j,:);
            Viterbi_Path_eq{i} = Viterbi_Path_eq{i}'; %%Gamma_Probs{i} = Gamma_Probs{i}';

        end

        % Transfer X_VEC to its 'real' values
        V_VEC_alpha_genotype = bitget(Viterbi_Path{i},1);
        V_VEC_beta_genotype = bitget(Viterbi_Path{i},5);
        V_VEC_alpha_copynumber = bitget(Viterbi_Path{i},2) + 2*bitget(Viterbi_Path{i},3) + 4*bitget(Viterbi_Path{i},4);
        V_VEC_beta_copynumber = bitget(Viterbi_Path{i},6) + 2*bitget(Viterbi_Path{i},7) + 4*bitget(Viterbi_Path{i},8);
        V_VEC_copy = V_VEC_alpha_copynumber+V_VEC_beta_copynumber;
        V_VEC_A = V_VEC_alpha_genotype.*V_VEC_alpha_copynumber + ...
            V_VEC_beta_genotype.*V_VEC_beta_copynumber;
        V_VEC_B = (1-V_VEC_alpha_genotype).*V_VEC_alpha_copynumber + ...
            (1-V_VEC_beta_genotype).*V_VEC_beta_copynumber;
        V_VEC_copy_eq = Viterbi_Path_eq{i};
        % Make a variable and force it to have copy numbers : 0 and 2
        V2 = Viterbi_Path{i}; V2(201:end) = V2(201:end)-2-32+64; 
        
            %  max_copy_alpha = max(V_VEC_alpha_copynumber)
        %  max_copy_beta = max(V_VEC_beta_copynumber)

        % % % % %         % Now plot the simulation results
        % % % % %         subplot(2,1,2);  hold on;
        % % % % %
        % % % % %         plot(V_VEC_alpha_genotype.*V_VEC_alpha_copynumber+0.06, 'g+');
        % % % % %         plot((1-V_VEC_alpha_genotype).*V_VEC_alpha_copynumber+0.06, 'go');
        % % % % %         plot(V_VEC_beta_genotype.*V_VEC_beta_copynumber-0.06, 'y+');
        % % % % %         plot((1-V_VEC_beta_genotype).*V_VEC_beta_copynumber-0.06, 'yo');
        % % % % %
        % % % % %         plot(V_VEC_alpha_copynumber-0.03, 'b');  plot(V_VEC_beta_copynumber+0.03, 'r');
        % % % % %         plot(V_VEC_alpha_copynumber+V_VEC_beta_copynumber+0.05, 'm');
        % % % % % %%       plot(V_VEC_A+0.2, 'ko'); plot(V_VEC_B+0.2, 'cx');
        % % % % % %%        plot(OUT_VEC, 'g.');  plot(OUT_VEC2, 'y.');
        % % % % %         legend('alpha copy number', 'beta copy number', 'SUM', 'A genotype', 'B genotpye', 'A Intensity', 'B Intensity');
        % % % % %         xlabel('SNP location'); ylabel('Values'); title('Viterbi reconstruction');


        % Now find the classification error
        if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
            CLASS_ERROR_MEAN(points_ind) = CLASS_ERROR_MEAN(points_ind) + sum( (V_VEC_copy - X_VEC_copy).^2 ./ length(X_VEC_copy) );
            CLASS_ERROR_STD(points_ind) = CLASS_ERROR_STD(points_ind) + (sum( (V_VEC_copy - X_VEC_copy).^2 / length(X_VEC_copy) ))^2;
            CLASS_ERROR_MEAN_eq(points_ind) = CLASS_ERROR_MEAN_eq(points_ind) + sum( (V_VEC_copy_eq - X_VEC_copy).^2 ./ length(X_VEC_copy) );
            CLASS_ERROR_STD_eq(points_ind) = CLASS_ERROR_STD_eq(points_ind) + (sum( (V_VEC_copy_eq - X_VEC_copy).^2 / length(X_VEC_copy) ))^2;

        else
            CLASS_ERROR_MEAN(points_ind) = CLASS_ERROR_MEAN(points_ind) + (sum(Viterbi_Path{i} ~= X_VEC) / length(X_VEC));
            CLASS_ERROR_STD(points_ind) = CLASS_ERROR_STD(points_ind) + (sum(Viterbi_Path{i} ~= X_VEC) / length(X_VEC))^2;
        end

        CLASS_ERROR_A_MEAN(points_ind) = CLASS_ERROR_A_MEAN(points_ind) + sum( (X_VEC_A-V_VEC_A).^2 ./ length(X_VEC_A) );
        CLASS_ERROR_B_MEAN(points_ind) = CLASS_ERROR_B_MEAN(points_ind) + sum( (X_VEC_B-V_VEC_B).^2 ./ length(X_VEC_B) );
        CLASS_ERROR_A_STD(points_ind) = CLASS_ERROR_A_STD(points_ind) + (sum( (X_VEC_A-V_VEC_A).^2 ./ length(X_VEC_A) ))^2;
        CLASS_ERROR_B_STD(points_ind) = CLASS_ERROR_B_STD(points_ind) + (sum( (X_VEC_B-V_VEC_B).^2 ./ length(X_VEC_B) ))^2;


        % Here the error is weighted according to the probability given
        BAYES_CLASS_ERROR_MEAN(points_ind) = BAYES_CLASS_ERROR_MEAN(points_ind) + (sum(abs(X_VEC'-Gamma_Probs{i}(2,:))) / length(X_VEC));
        BAYES_CLASS_ERROR_STD(points_ind) = BAYES_CLASS_ERROR_STD(points_ind) + (sum(abs(X_VEC'-Gamma_Probs{i}(2,:))) / length(X_VEC))^2;

        % Here the error is based on 0.5 threshold clipping !!!
        %       Gamma_Clips =  (Gamma_Probs{i}(2,:) > 0.5);
        %       BAYES_CLASS_ERROR_MEAN(points_ind) = BAYES_CLASS_ERROR_MEAN(points_ind) + (sum(abs(X_VEC'-Gamma_Clips)) / length(X_VEC));
        %       BAYES_CLASS_ERROR_STD(points_ind) = BAYES_CLASS_ERROR_STD(points_ind) + (sum(abs(X_VEC'-Gamma_Clips)) / length(X_VEC))^2;

        figure; subplot(2,2,1); plot(X_VEC_alpha_copynumber, '.'); legend('\alpha copy');    
        subplot(2,2,2); plot(X_VEC_beta_copynumber, '.'); legend('\beta copy');            
        subplot(2,2,3); plot(V_VEC_alpha_copynumber, '.'); legend('recovered \alpha copy');    
        subplot(2,2,4); plot(V_VEC_beta_copynumber, '.'); legend('recovered \beta copy');                  
    end % iter loop

    points_ind = points_ind + 1;
    done_num_points = num_points
end  % num_points loop

% Normalize the errors
M_SQR_ERROR_MEAN = M_SQR_ERROR_MEAN ./ num_iters;
MEW_SQR_ERROR_MEAN = MEW_SQR_ERROR_MEAN ./ num_iters;
SIGMA_SQR_ERROR_MEAN = SIGMA_SQR_ERROR_MEAN ./ num_iters;
KL_ERROR_MEAN = KL_ERROR_MEAN ./ num_iters;
CLASS_ERROR_MEAN = CLASS_ERROR_MEAN ./ num_iters;
BAYES_CLASS_ERROR_MEAN = BAYES_CLASS_ERROR_MEAN ./ num_iters;
CLASS_ERROR_A_MEAN = CLASS_ERROR_A_MEAN ./ num_iters;
CLASS_ERROR_B_MEAN = CLASS_ERROR_B_MEAN ./ num_iters;
CLASS_ERROR_MEAN_eq = CLASS_ERROR_MEAN_eq ./ num_iters;



M_SQR_ERROR_STD = sqrt(M_SQR_ERROR_STD ./ num_iters - (M_SQR_ERROR_MEAN .^ 2));
MEW_SQR_ERROR_STD = sqrt(MEW_SQR_ERROR_STD ./ num_iters - (MEW_SQR_ERROR_MEAN .^ 2));
SIGMA_SQR_ERROR_STD = sqrt(SIGMA_SQR_ERROR_STD ./ num_iters - (SIGMA_SQR_ERROR_MEAN .^ 2));
KL_ERROR_STD = sqrt(KL_ERROR_STD ./ num_iters - (KL_ERROR_MEAN .^ 2));
CLASS_ERROR_STD = sqrt(CLASS_ERROR_STD ./ num_iters - (CLASS_ERROR_MEAN .^ 2));
BAYES_CLASS_ERROR_STD = sqrt(BAYES_CLASS_ERROR_STD ./ num_iters - (BAYES_CLASS_ERROR_MEAN .^ 2));
CLASS_ERROR_A_STD = sqrt(CLASS_ERROR_A_STD ./ num_iters - (CLASS_ERROR_A_MEAN .^ 2));
CLASS_ERROR_B_STD = sqrt(CLASS_ERROR_B_STD ./ num_iters - (CLASS_ERROR_B_MEAN .^ 2));
CLASS_ERROR_STD_eq = sqrt(CLASS_ERROR_STD_eq ./ num_iters - (CLASS_ERROR_MEAN_eq .^ 2));

% Now do the plots

colorvec = 'bgrkmc';

% Show learned model parameters
if(learn_flag)
    figure; hold on;  title('Relative Errors in M - Transition Matrix'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    % Plot M errors
    for i=1:HMM_x_dim
        for j=1:HMM_x_dim
            errorbar(num_points_vec, reshape(M_SQR_ERROR_MEAN(i,j,:), 1, length(num_points_vec)), reshape(M_SQR_ERROR_STD(i,j,:), 1, length(num_points_vec)), colorvec(i*HMM_x_dim+j));
        end
    end

    figure; hold on; title('Relative Errors in MEW (divided by mew gaps) - Expectation Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    % Plot MEW errors
    for i=1:HMM_x_dim
        errorbar(num_points_vec, MEW_SQR_ERROR_MEAN(i,:), MEW_SQR_ERROR_STD(i,:), colorvec(i));
    end

    figure; hold on; title('Relative Errors in SIGMA - Standard Error Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error');
    % Plot MEW errors
    for i=1:HMM_x_dim
        errorbar(num_points_vec, SIGMA_SQR_ERROR_MEAN(i,:), SIGMA_SQR_ERROR_STD(i,:), colorvec(i));
    end


    % Plot the KL distance
    figure; hold on; title('Kullback-Leibler Error Vector'); xlabel('Num. Samples'); ylabel('KL Error');
    % size(num_points_vec)
    % size(KL_ERROR_MEAN)
    % size(KL_ERROR_STD)
    errorbar(num_points_vec, KL_ERROR_MEAN, KL_ERROR_STD);
end


if(HMM_MODEL.SPECIAL_MODELS_FLAG == 1)
    % % % %     % Plot the four errors
    % % % %     figure; subplot(2,2,1); hold on; title('alpha A');
    % % % %     plot( (X_VEC_alpha_genotype.*X_VEC_alpha_copynumber - V_VEC_alpha_genotype.*V_VEC_alpha_copynumber).^2, '.');
    % % % %     subplot(2,2,2); hold on; title('alpha B');
    % % % %     plot( ((1-X_VEC_alpha_genotype).*X_VEC_alpha_copynumber - (1-V_VEC_alpha_genotype).*V_VEC_alpha_copynumber).^2, '.');
    % % % %     subplot(2,2,3); hold on; title('beta A');
    % % % %     plot( (X_VEC_beta_genotype.*X_VEC_beta_copynumber - V_VEC_beta_genotype.*V_VEC_beta_copynumber).^2, '.');
    % % % %     subplot(2,2,4); hold on; title('beta B');
    % % % %     plot( ((1-X_VEC_beta_genotype).*X_VEC_beta_copynumber - (1-V_VEC_beta_genotype).*V_VEC_beta_copynumber).^2, '.');

    figure; subplot(2,1,1); hold on; plot( (X_VEC_A - V_VEC_A).^2-0.05, '.');
    plot( (X_VEC_B - V_VEC_B).^2+0.05, 'r.'); plot( (V_VEC_copy-X_VEC_copy).^2+0.1,'g.'); 
    plot( (V_VEC_copy_eq-X_VEC_copy).^2, 'm.');
    legend('A error', 'B error', 'Tot. error', 'Tot. error eq.');
    subplot(2,1,2); hold on; plot( (X_VEC_alpha_genotype - V_VEC_alpha_genotype).^2, '.');
    plot( (X_VEC_beta_genotype - V_VEC_beta_genotype).^2, 'r.');
    legend('\alpha error', '\beta error');

    figure; hold on; title('Classification (Viterbi) Error Probability Vector A,B');
    xlabel('Num. Samples'); ylabel('CLASS. Error');
    errorbar(num_points_vec, CLASS_ERROR_A_MEAN, CLASS_ERROR_A_STD);
    errorbar(num_points_vec+1, CLASS_ERROR_B_MEAN, CLASS_ERROR_B_STD,'r');
    errorbar(num_points_vec+2, CLASS_ERROR_MEAN, CLASS_ERROR_STD,'g');
    errorbar(num_points_vec+3, CLASS_ERROR_MEAN_eq, CLASS_ERROR_STD_eq,'m');
    legend('A', 'B', 'Total', 'Total eq.');

else
    % Plot the Classification (Viterbi) error
    figure; subplot(2,1,1); hold on; title('Classification (Viterbi) Error Probability Vector'); xlabel('Num. Samples'); ylabel('CLASS. Error');
    errorbar(num_points_vec, CLASS_ERROR_MEAN, CLASS_ERROR_STD);
    % Plot the Bayesian Classification (forward) error
    subplot(2,1,2); hold on; title('Classification (forward) "Error Probability" Vector'); xlabel('Num. Samples'); ylabel('forward CLASS. Error');
    errorbar(num_points_vec, BAYES_CLASS_ERROR_MEAN, BAYES_CLASS_ERROR_STD);
end


% Here if we want we plot the real data and the Viterbi chosen path
% % % figure; hold on;
% % %
% % % plot( OUT_VEC(1:200), '*m');
% % % plot( Viterbi_Path{i}(1:200), 'bx');
% % % ylabel('Amp. Level ', 'fontsize', 8);   % Do the plot
