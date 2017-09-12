% Here witerse load a data, computec the correlations and geOt the desired
% fraction
path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
path(path, 'E:\Research\PACPerm\numeric');
% path(path, 'C:\Weizmann\Research\PACPerm\numeric');
% path(path, 'C:\Weizmann\Research\PACPerm\numeric\nmm\linalg');
%
% path(path, '/a/fs-01/pdoc/mrosenz/Liat/GeneCorr/Or');
path(path, 'E:\Gene_Corr\');

ttt_time = cputime;
TRUE = 1; FALSE = 0;

TOL = 0.00000001;

kept_fraction=0.01;
KNN = 10; % parameter for knn missing values algorithm
Fisher_R_to_Z = 1; % Flag saying if to perform fisher RtoZ transformation
remove_zero_flag = TRUE;   % Flag saying to remove all the features with correlation
% zero from the data, since they cause bias.


check_eps_delta_flag = 1;

% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; FROM_DATA = 3; mix_GAUSSIAN=4;
student_t=5;% in the last one we take q simply according to the data bins
rand_flag = mix_GAUSSIAN;

% data files
OLD_VANT_VEER = 1; NEW_ROSSETA = 2; WANG = 3;     % gene expression
TOPIC = 4; % Author-topic matching
NIPS_DOROTHEA = 5; NIPS_ARCENE = 6; NIPS_DEXTER = 7; NIPS_GISETTE = 8; % Various datasetes from NIPS contest
ROSEN = 9; HEPATO =10; % New data's from Assif
Bhattacharjee=11;  YEOH = 12;  % New Lung data
BRAIN = 12; KIM = 13; % New aging related datasets

RAND_DATA = 14; % Here we randomize a data so that we have control on it, and also can work with matlab 6.5 (no loading from files)


% Vector saying which datasets come in a sparse form
IS_SPARSE_VEC =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0];
CALC_CORRS_VEC = [1,1,1,0,0,0,0,1,1,1,1,1,1,1];  % calc for Gizette [1,1,1,0,0,0,0,0];
%%%CALC_VAR_VEC = [0,0,0,0,0,0,0,0,0];
CALC_VAR_VEC =   [1,1,1,1,1,1,1,1,1,1,1,1,1,1];

calc_corrs_flag = FALSE; % flag saying if to calculate the correlations from the data, or do something else

% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TWO_SAMPLED;



%alpha = 0.014;
alpha_vec =  [ 0.12 0.00403001 0.012]; % corresponds to ~70 and ~700 genes
%TOP_GENES = 70;

epsilon_needed = 0.7; delta_needed = 0.5; % Desired error and confidence rates


alpha=alpha_vec(1);
miu=[];
prior=[];
%nsamples=3;
iters = 100;



% Now chose which data to load - New : Run on all data's !!!!

num_samples_needed = zeros(1,2);
std_needed = zeros(1,2);


SMALL_GENES_NUM = 1000; % Take only this amount of genes ...

for data_flag =  HEPATO % YEOH % HEPATO  % NEW_ROSSETA; % BRAIN; % BRAIN; % HEPATO % OLD_VANT_VEER % HEPATO % TOPIC % OLD_VANT_VEER % NEW_ROSSETA % NIPS_GISETTE % NEW_ROSSETA % NIPS_GISETTE; %NIPS_DEXTER; %%NIPS_ARCENE; %NIPS_DOROTHEA; % TOPIC; % WANG; %%% NEW_ROSSETA; %%%OLD_VANT_VEER;
    cur_data_is = data_flag
    Fisher_R_to_Z = 1; % Flag saying if to perform fisher R-to-Z transformation
    one_side_flag = TWO_SIDES; % make sure we do two-sided unless later we change the flag
    calc_corrs_flag = CALC_CORRS_VEC(data_flag);
    estimate_variance_flag = CALC_VAR_VEC(data_flag);
    % Now for each data we do one figure !
    plot_ind=1;
    %%%%    figure; %subplot(2,2,1);     hold on;

    clear R; % free memory


    if(data_flag == OLD_VANT_VEER)
        rand_flag = GAUSSIAN; % should be gaussian
        num_of_Gaussians = 2;
        load('../data/breastData.mat'); % The name of the data file containing the expression and labels.
        data_str = 'VantVeerOld';
    end
    if(data_flag == NEW_ROSSETA)
        loading_big_rosseta = 1;
        load('../data/Rosetta_data.mat'); % The name of the data file containing the expression and labels.
        data_str = 'VantVeerNew';

        finished_loading_big_rosseta = 1
        one_side_flag = TWO_SIDES;
        remove_zero_flag = FALSE;
        rand_flag = GAUSSIAN;
        %%%%%%%%removing genes with too many missing values
        present=real(pval<0.05);
        present_idx=find(sum(present,2)>=size(pval,2)*kept_fraction);
        present_data=data(present_idx,:);
        data=present_data;
        R.dat = knnimpute(data, KNN); %         R.dat = data;
        calc_corrs_flag = TRUE;
        %         R.dat = data;
        rand_flag = GAUSSIAN;
        R.Labels = real(Labels<100); R.Labels = R.Labels';
        loading_big_rosseta = 1
    end

    if(data_flag == WANG)
        load('../data/WANG_DATA.mat');
        data_str = 'Wang';

        rand_flag = GAUSSIAN;
        %%%%%%%removing genes with too many missing values
        present=real(pval<0.05);
        present_idx=find(sum(present,2)>=size(pval,2)*kept_fraction);
        present_data=gene_expression_log2(present_idx,:);
        data=present_data;
        % Here we must do log since Gary didn't do it
        R.dat = data;
        R.Labels = Prognosis';
        calc_corrs_flag = TRUE;
    end

    if(data_flag == TOPIC)
        load('../data/topic_data.mat');
        data_str = 'Author-Topic';
        rand_flag = mix_GAUSSIAN;
        R.dat = Topics_N_counts;
        R.Labels = N_Labels';
        Pearson_corrs = words_corr;   % Many zero correlations in this dataset.

        % avoid constant rows by randomizing them
        const_rows_indexes = find(max(R.dat,[],2) == min(R.dat,[],2));
        %        R.dat(const_rows_indexes,:) = rand(length(const_rows_indexes), size(R.dat,2));
        num_of_Gaussians = 2;

    end

    if(data_flag == NIPS_DOROTHEA)
        load('../data/an_dorothea.mat');
        data_str = 'NIPS-Dorothea';
        R.dat = knnimpute(T, KNN);   %%%%%   R.dat = T;
        R.Labels = Labels';
        Pearson_corrs = corr;   % Note : Here these aren't correlations but something else ! Distribution not gaussian but a lot of mass on the right side
        one_side_flag = ONE_SIDE; % here we take only positive correlations
        Fisher_R_to_Z = 0; % Do not do Fisher's transform since we do not deal with corrleations here
    end

    if(data_flag == NIPS_ARCENE)
        load('../data/ARCENE_data_file.mat');
        data_str = 'NIPS-Arcene';
        R.dat = Training.data;
        R.Labels = Labels';
        Pearson_corrs = corr;   % Very weak correlations in this dataset
    end

    if(data_flag == NIPS_DEXTER)
        load('../data/DEXTER_data_file.mat');
        data_str = 'NIPS-Dexter';
        R.dat = T;
        R.Labels = Labels';
        Pearson_corrs = corr;   % Very weak correlations in this dataset. Most are zero. NOT a gaussian distribution
    end

    if(data_flag == NIPS_GISETTE)
        load('../data/GISETTE_data_file.mat');
        data_str = 'NIPS-Gisette';
        rand_flag = GAUSSIAN; %%% mix_GAUSSIAN; not supported yet ..
        num_of_Gaussians=2;
        num_of_itterations=100;

        %         rand_flag =student_t;
        R.dat = T;

        % avoid constant rows by randomizing them
        const_rows_indexes = find(max(R.dat,[],2) == min(R.dat,[],2));
        R.dat(const_rows_indexes,:) = rand(length(const_rows_indexes), size(R.dat,2));

        %         start_loop = 2
        %         % Correct R.dat to not have any non-numeric
        %         for i=1:size(R.dat,1)
        %             for j=1:size(R.dat,2)
        %                 if(isnumeric(R.dat(i,j)==0))
        %                     R.dat(i,j) = rand(1);
        %                 end
        %             end
        %         end
        %
        %         end_loop = 3

        R.Labels = Labels';
        Pearson_corrs = corr;   % Very weak correlations in this dataset. Not symmetric, mostly negative correlations.
    end

    if(data_flag == ROSEN)
        load('../data/ROSEN.mat');
        data_str = 'ROSEN';
        rand_flag = GAUSSIAN;
        R.dat = data;
        R.Labels = labels;   % A gaussian distribution of correlations (?)

        % Now take only a small subset of the genes !
        %temp_Gary = sum(R.dat < 0.5,2);
        %temp_Ind = find(temp_Gary < 200);

        %R.dat = R.dat(temp_Ind,:);

        % do knn impute, threshold to 0.5 and take log2 !!!
        R.dat = knnimpute(R.dat, KNN);

        % don't do log - this is already after log!
        % R.dat = log2(R.dat);

        %      R.dat = R.dat(1:SMALL_GENES_NUM,:);
    end

    if(data_flag == HEPATO)
        load('../data/hepatocellular_carcinoma.mat');
        data_str = 'Hepatocellular-Carcinoma';
        rand_flag = GAUSSIAN;
        R.dat = data;
        R.Labels = labels'; % A gaussian distribution of correlations (?)

        R.dat(R.dat < 20) = 20;
        R.dat = R.dat + 0.001*rand(size(R.dat,1),size(R.dat,2)); % Add very small random noise
        
        temp_Gary = sum(R.dat < 30,2);
        temp_Ind = find(temp_Gary < 15);

        R.dat = R.dat(temp_Ind,:);
        R.dat = knnimpute(R.dat, KNN); % This is now done over negative numbers which is baddd.
        R.dat = log2(R.dat);

    end

    if(data_flag == BRAIN)
        load('../data/brain_age.mat');
        data_str = 'Aging-Brain';
        rand_flag = GAUSSIAN;

        %   R.dat(R.dat < 30) = 30;
        temp_Gary = sum(R.dat < 30,2);
        temp_Ind = find(temp_Gary < 15);

        R.dat = R.dat(temp_Ind,:);
        R.dat = knnimpute(R.dat, KNN);
        R.dat = log2(R.dat);

        % R.dat = R.dat(1:SMALL_GENES_NUM,:);

    end
    
    

   if(data_flag == Bhattacharjee)
        rand_flag = GAUSSIAN;corr_
        load('../data/Bhattacharjee_66.mat'); % The name of the data file containing the expression and labels.
        data_str = 'Bhattacharjee';
        data(data<30)=30;
        data = data + 0.001*randn(size(data,1),size(data,2));
        data=log2(data);
        Y = var(data,0,2);
        [val,ind]=sort(-Y);
        present_data=data(ind(1:4000),:); % Take highest variance genes
        R.dat=present_data;
        R.Labels=labels; 
   end

   
   
   if(data_flag == YEOH)
       rand_flag = GAUSSIAN;
       kept_fraction=.9;
       d_var_n=0;
       load('../data/yeoh'); % The name of the data file containing the expression and labels.
       res = 4; % the resulution of N that we use
       max_nsamples = floor(size(data,2)/2); % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
       min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
       nsamples_vec = min_nsamples:res:max_nsamples;
       data_str = 'YEOH';
       pval=(data<30);
       data(data<30)=30;
       data=log2(data);
       present=real(pval<0.05);
       N_s=size(data,2);

       present_idx=find(sum(present,2)>=N_s*kept_fraction);
       present_data=data(present_idx,:);
       R.dat=present_data;
       R.Labels=labels;
   end


   
   if(data_flag == KIM)
       load('../data/KimKidney.mat');
       data_str = 'Aging-Kidney';
       rand_flag = GAUSSIAN;

       %   R.dat(R.dat < 30) = 30;
       temp_Gary = sum(R.dat < 1,2);
       temp_Ind = find(temp_Gary < 15);

       R.dat = R.dat(temp_Ind,:);

       R.dat(R.dat<1)=1;
       R.dat = R.dat + 0.001*randn(size(R.dat,1),size(R.dat,2));

       R.dat = log2(R.dat);

       % R.dat = R.dat(1:SMALL_GENES_NUM,:);

   end


    if(data_flag == RAND_DATA)
        
        DISCRETE_LABELS = 1; CONTINUOUS_LABELS = 0;% 1 discrete 0/1 0 Gaussians

        labels_type = DISCRETE_LABELS;

        Ngenes = 5000; Nsamples = 30;

        if(labels_type == DISCRETE_LABELS)
            R.Labels = round(rand(1, Nsamples)); % 0/1 labels
            pos_ind = find(R.Labels == 1);
        else
            R.Labels = randn(1, Nsamples); % 0/1 labels
        end

        rand_flag = GAUSSIAN;


        % Generate a more sophisticated  random model: % Here small
        % correlations mean no signal!
        Q_corr_vec = 0.02*(randn(1, Ngenes)); % Correlation of each gene with the survival
        W_corr_vec = 0.02*(randn(1, Ngenes));

        % Now for each gene we need to generate an expression value with the given correlation. Maybe by MOG?
        if(labels_type == DISCRETE_LABELS) % MOG
            mu_vec = 2.*tanh(Q_corr_vec) ./ sqrt(1 - tanh(Q_corr_vec).^2);

            R.dat = randn(Ngenes,Nsamples);
            R.dat(:,find(R.Labels==1)) = R.dat(:,find(R.Labels==1)) + repmat(mu_vec, length(pos_ind), 1)';
        else
            R.dat = randn(Ngenes,Nsamples); % First generate the r.v.s

            C_corr_vec = tanh(Q_corr_vec);
            beta_vec = sqrt(C_corr_vec .^ 2 ./ (1 - C_corr_vec .^ 2)) .* sign(C_corr_vec);

            % Add to generate the correlations:
            R.dat = R.dat + (repmat(R.Labels, Ngenes, 1) .* repmat(beta_vec', 1, Nsamples));

        end % if part

        data_str = 'RAND-data';
    end


    Ngenes = size(R.dat, 1)
    Nsamples = size(R.dat, 2)


    if(calc_corrs_flag == TRUE)
        % Calculate the mean correlation of each gene with survival
        normed_labels = R.Labels - mean(R.Labels); normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
        normed_data = R.dat;
        %        normed_data(4965,:) = rand(1,Nsamples); % BADD THING !!!
        normed_data = normed_data - repmat(mean(normed_data, 2), 1, Nsamples);
        normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, Nsamples);
        %normed_data = normed_data';

        Pearson_corrs = sum(normed_data .* repmat(normed_labels, Ngenes, 1),2);
        size(Pearson_corrs)

        % Generate the correlation matrix
        corr_mat = [];

        % % % % % % %         corr_mat = normed_data * normed_data';
        % % % % % % %
        % % % % % % %         % Sor the correlation matrix of the genes according to their correlation with outcome
        % % % % % % %         [sorted_vec sort_inds ] = sort(Pearson_corrs);
        % % % % % % %         corr_mat = corr_mat(sort_inds,sort_inds);


    else % here the correlations are already given

        if(remove_zero_flag == TRUE)  % remove zero correlations only if we do not calculate them
            non_zero_ind = find(Pearson_corrs ~= 0);
            Pearson_corrs = Pearson_corrs(non_zero_ind);
            R.dat = R.dat(non_zero_ind,:);
            Ngenes = size(R.dat, 1)


            % avoid 'almost' constant rows by randomizing them
            const_rows_indexes = find(max(R.dat,[],2) == min(R.dat,[],2));

            R.dat = R.dat + TOL * rand(size(R.dat,1), size(R.dat,2));

        end
    end
    index = 1;


    max_nsamples = Nsamples*2; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
    small_res = floor(Nsamples/30); res = 4*small_res;  max_nsamples = 15*res;   % Set the resolution we want!
    nsamples_vec = [  res:small_res:15*small_res]; small_vec_len = 15-3 %  3*res:res:max_nsamples]; small_vec_len = 15-3;
    res
    small_res

    % Choose if you want to have short vectors !
    nsamples_vec = nsamples_vec(9:end); small_vec_len = length(nsamples_vec); % Semi-Short Vecs 
   %%%% nsamples_vec=nsamples_vec(end); small_vec_len = 1; % Short vecs !!! 

    %     nsamples_vec = [  small_res*8:small_res:20*small_res ]; small_vec_len = length(nsamples_vec);

    % % % % %     % Now do the plot that shows if indeed the Fisher Z's have std which is
    % % % % %     % independent of their value
    % % % % %     if( (calc_corrs_flag == TRUE) & (estimate_variance_flag == TRUE))
    % % % % %         corrs_iters = 2000;
    % % % % %         [N_VEC Pearson_corrs Fisher_Zs Rho_sig_mean Rho_sig_std Rho_bias_mean Rho_bias_std Z_sig_mean Z_sig_std Z_bias_mean ...
    % % % % %             Z_bias_std chunk_Rho_sig chunk_Z_sig chunk_Rho_bias chunk_Z_bias ] = ...
    % % % % %             CalcZVarianceFromDataFunc(R, alpha_vec, rand_flag, corrs_iters ); % take only max samples to save time
    % % % % %         figure; %subplot(2,3,4);
    % % % % %         hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. N=' num2str(N_VEC)]);
    % % % % %         xlabel('true P'); ylabel('std');
    % % % % %         figure; %subplot(2,3,5);
    % % % % %         hold on; plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. N='  num2str(N_VEC(1))]);
    % % % % %         xlabel('true Z'); ylabel('std');
    % % % % %     end   % plot of std of the



    % Now call the function to do all work on estimating the variance and
    % distribution of the Z scores
    if(  (estimate_variance_flag == TRUE))  %%% (calc_corrs_flag == TRUE) &  % For now we do this always - not demand that we calculate corrs
        corrs_iters = 10;
        [N_VEC Pearson_corrs Fisher_Zs Rho_sig_mean Rho_sig_std Rho_bias_mean Rho_bias_std Z_sig_mean ...
            Z_sig_std Z_bias_mean Z_bias_std chunk_Rho_sig chunk_Z_sig chunk_Rho_bias chunk_Z_bias hist_Pearson_Ps hist_Fisher_Zs random_genes_picked] = ...
            CalcZVarianceFromDataFunc(R, alpha, rand_flag, corrs_iters);
        
        corrs_iters2=10;

        
        % Skip the following: Very time-consuming
        [N_VEC2 N_HALF_VEC2 Z_sig_mean_plot Q_sig_mean_plot ] = CalcTwoZVariancesFromDataFunc(R, alpha, rand_flag, corrs_iters2 );
       
        
         % Plot and compare both variances
         figure; hold on; plot(N_HALF_VEC2, Z_sig_mean_plot, '.'); plot(N_VEC2, Q_sig_mean_plot, '.r');  plot(N_HALF_VEC2, Q_sig_mean_plot(1:length(N_HALF_VEC2))-Z_sig_mean_plot, '.g'); 
         plot(N_HALF_VEC2, 1 ./ sqrt(N_HALF_VEC2-3), 'k');    plot(N_VEC2, std(Fisher_Zs) + 1 ./ sqrt(N_VEC2-3), 'k');
         legend('Ave Z std.', 'Q std.', 'diff', 'Fisher', 'Fisher plus one'); xlabel('Nsamples'); ylabel('Variance'); 
         title(['Std. of each genes Z and of total Q, Data is ' data_str]);
    end

    
    % Make some dummy variables
    if(data_flag ~= RAND_DATA)
        Q_corr_vec = Fisher_Zs;
        W_corr_vec = Fisher_Zs;
    end

    % Try to fit the std. of Z according to its mean :
    quad_fit =  fittype('a1*x^2+a2*x+a3');
    %%%%%    POLY = fit(Fisher_Zs, chunk_Z_sig', quad_fit);
    POLY = fit(Fisher_Zs, 2*mean(chunk_Z_sig')-chunk_Z_sig', quad_fit);
    coeffs_vec = [POLY.a1 POLY.a2 POLY.a3] .* sqrt(Nsamples/2); % Adjust for the number of samples !


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot to see correlation between the 'true' value and the variance
    min_x_val = min(Fisher_Zs) - 0.01; max_x_val = max(Fisher_Zs) + 0.01;
    min_y_val = min( min(chunk_Z_sig), min(chunk_Rho_sig) ) - 0.01;
    max_y_val = max( max(chunk_Z_sig), max(chunk_Rho_sig) ) + 0.01;
    figure; subplot(2,3,1); hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. n_{samp}=' num2str(N_VEC(end))]);
    xlabel('true P'); ylabel('std'); axis([min_x_val max_x_val min_y_val max_y_val]);
    subplot(2,3,4); hold on; plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. n_{samp}='  num2str(N_VEC(end))]);
    xlabel('true Z'); ylabel('std'); axis([min_x_val max_x_val min_y_val max_y_val]);
    min_y_val = min( min(chunk_Z_bias), min(chunk_Rho_bias) ) - 0.01;
    max_y_val = max( max(chunk_Z_bias), max(chunk_Rho_bias) ) + 0.01;
    subplot(2,3,2); hold on; plot(Pearson_corrs, chunk_Rho_bias, '.'); title(['True Pearson P and estimator bias.  n_{samp}=' num2str(N_VEC(end))]);
    xlabel('true P'); ylabel('bias'); axis([min_x_val max_x_val min_y_val max_y_val]);
    subplot(2,3,5); hold on; plot(Fisher_Zs, chunk_Z_bias, '.'); title(['True Fisher Z and estimator bias. n_{samp}='  num2str(N_VEC(end))]);
    xlabel(['true Z, for n_{samp} = ' num2str(N_VEC(end)), ' for data: ' data_str]); ylabel('bias'); axis([min_x_val max_x_val min_y_val max_y_val]);

    % Plot to see the histogram of the variance of Z
    subplot(2,3,3); hold on; hist(chunk_Rho_sig, 100); title(['True Pearson P std. hist. n_{samp}=' num2str(N_VEC(end))]);
    xlabel('std. P'); ylabel('Freq. ');
    subplot(2,3,6); hold on; hist(chunk_Z_sig, 100); title(['True Fisher Z std. hist. n_{samp}='  num2str(N_VEC(end))]);
    xlabel('std. Z'); ylabel('Freq. ');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot for a subset of 16 randomly selected genes, the histogram distribution of the Fisher Z score on a finite sample of size ???
    figure; subplot(4,4,1); hold on;
    min_x_val = min(min(hist_Fisher_Zs)) - 0.01;
    max_x_val = max(max(hist_Fisher_Zs)) + 0.01;
    for i=1:4
        for j=1:4
            subplot(4,4,4*(i-1)+j); hold on;
            h = hist(hist_Fisher_Zs(4*(i-1)+j,:), 100);
            hist(hist_Fisher_Zs(4*(i-1)+j,:), 100); title(['gene # ' num2str(random_genes_picked(4*(i-1)+j))]);
            axis([min_x_val max_x_val 0 max(h)+1] );
        end
    end
    subplot(4,4,15); hold on; xlabel(['Fisher Zs for n_{samp} = ' num2str(N_VEC(end)) ' , ' ] );
    subplot(4,4,16); hold on; xlabel([' Data is: ' data_str ] );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot for a subset of 16 randomly selected genes, the histogram distribution of the Fisher Z score on a finite sample of size ???
    figure; subplot(4,4,1); hold on;
    min_x_val = min(min(R.dat)) - 0.01;
    max_x_val = max(max(R.dat)) + 0.01;
    for i=1:4
        for j=1:4
            subplot(4,4,4*(i-1)+j); hold on;
            h = hist(R.dat(random_genes_picked(4*(i-1)+j),:), 100);
            hist(R.dat(random_genes_picked(4*(i-1)+j),:), 100);     title(['gene # ' num2str(random_genes_picked(4*(i-1)+j))]);
            axis([min_x_val max_x_val 0 max(h)+1] );
        end
    end
    subplot(4,4,15); hold on; xlabel(['Values of various features (genes) , '] );
    subplot(4,4,16); hold on; xlabel([' Data is: ' data_str ] );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    % Remove the mean
    normed_hist_Fisher_Zs = hist_Fisher_Zs - repmat(Fisher_Zs(random_genes_picked), 1, corrs_iters);

    % Choose 16 random couples
    figure; subplot(4,4,1); hold on; title('See noises correlations');
    rand_second_ind = randperm(16);
    for i=1:16
        subplot(4,4,i); hold on;
        plot(normed_hist_Fisher_Zs(i,:),normed_hist_Fisher_Zs(rand_second_ind(i),:),'*');
    end


    % Calculate all correlations and make an histogram
    NUM_RAND_GENES = 36;



    % Calculate also the original correlation matrix
    Orig_CorrMatrix = zeros(NUM_RAND_GENES);
    for i=1:NUM_RAND_GENES
        for j=1:NUM_RAND_GENES
            Orig_CorrMatrix(i,j) = corr(R.dat(random_genes_picked(i),:)', R.dat(random_genes_picked(j),:)');
        end
    end

    % Calculate the correlation to the Survival vector and display !!!
    for i=1:NUM_RAND_GENES
        Orig_CorrMatrix(i,end) = corr(R.dat(random_genes_picked(i),:)', R.Labels');
        Orig_CorrMatrix(end,i) = Orig_CorrMatrix(i,end);
    end
    Orig_CorrMatrix(end,end) = 1;

    [sorted_vec sort_inds ] = sort(Orig_CorrMatrix(:,end));



    figure; imagesc(Orig_CorrMatrix(sort_inds,sort_inds)); colorbar; title('Correlations in Expression');


    % Matrix of the correlations in Z's noise
    Var_Z_CorrMatrix = zeros(NUM_RAND_GENES);

    for i=1:NUM_RAND_GENES
        for j=1:NUM_RAND_GENES
            Var_Z_CorrMatrix(i,j) = corr(normed_hist_Fisher_Zs(i,:)', normed_hist_Fisher_Zs(j,:)');
        end
    end


    figure; imagesc(Var_Z_CorrMatrix(sort_inds(1:end-1),sort_inds(1:end-1))); colorbar; title('Correlations in noise of correlation with survival');



    figure; title('Histogram of correlation of genes Zs noise'); hist(Var_Z_CorrMatrix(:), 100);

    figure; title('Histogram of correlation of genes Expression'); hist(Orig_CorrMatrix(:), 100);



    figure; plot(Orig_CorrMatrix(:), Var_Z_CorrMatrix(:), '.');
    hold on; xlabel('orig exp.'); ylabel('noise expr'); title('Expression correlation vs. Z noise correlation');


   %%%return;

   % Run on different iterations
   iters_vec = [100,500,1000,2000];
   
for iters = iters_vec    
    now_iter_is = iters
    for true_corr_flag=0:0 % Go over the two options   %%% TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;

        % Now try to estimate the variance of the correlation of each gene with
        % survival. The method we use is go over all the couples.
        % Below the 2-factor is for both values of true_corr_flag
        inf_limit_frac = zeros(2 * length(alpha_vec), length(nsamples_vec));
        inf_limit_std = zeros(2 * length(alpha_vec), length(nsamples_vec));
        samp_frac_mean = zeros(2 * length(alpha_vec), length(nsamples_vec));
        samp_frac_std = zeros(2 * length(alpha_vec), length(nsamples_vec));
        samp_frac_mean_from_data = zeros(2 * length(alpha_vec), length(nsamples_vec));
        samp_frac_std_from_data = zeros(2 * length(alpha_vec), length(nsamples_vec));

        %%%     for beta = alpha_vec % alpha_vec alpha
        alpha = alpha_vec(1); % TOP_GENES / Ngenes; % make the number of genes always constant (70)
        true_flag_is = true_corr_flag
        alpha_is = alpha
        start_plot_ind = plot_ind

        if(Fisher_R_to_Z)

            % Do the hyperbolic tangent transform :
            Fisher_Zs = 0.5 * (log(1+Pearson_corrs) - log(1-Pearson_corrs));
        else
            % Here the names are misleading and do not descibe the content
            % of the variables :
            Fisher_Zs = Pearson_corrs;
        end

        SAVE_FISHER_Z{data_flag} = Fisher_Zs; % Save the Fisher Z's for plotting later - Note: This is calculated for all the data (Nsamples)

        % Now try to fit the best linear/gaussian/whatever distribution to the correlations data.
        % We try with the three options to fit : gaussian, uniform and linear.
        % So far only Gaussian is supported
        if(rand_flag == GAUSSIAN)
            Fisher_Zs_mean = mean(Fisher_Zs);
            Fisher_Zs_std = std(Fisher_Zs);
            SAVE_FISHER_Z_MEAN{data_flag} = Fisher_Zs_mean;
            SAVE_FISHER_Z_STD{data_flag} = Fisher_Zs_std;

            % Dummy for MOG
            miu = Fisher_Zs_mean; prior = 1;

            Fisher_Zs_std
            % New ! Make a correction to the std. in order to
            % compensate for the noise of the samples
 %%%%           Fisher_Zs_std = sqrt(Fisher_Zs_std^2 -  (1/(Nsamples-1*3)))  % TRY FIT - TEMPORARY !!!! Alternative : Estimate from data!
  %%%%          Fisher_Zs_std = sqrt(Fisher_Zs_std^2 -
  %%%%          Z_sig_mean(end)^2/2);  % Here we take the sigma from the
  %%%%          data and don't trust the 1/(n-3) !!!
            Fisher_Zs_std = 1; % Temp, for RAND-DATA we know the variance of Q
  
            Fisher_Zs_min = min(Fisher_Zs);  Fisher_Zs_max = max(Fisher_Zs); Fisher_Zs_gap = Fisher_Zs_max-Fisher_Zs_min;

            Fisher_Z_size = size(Fisher_Zs)

            Fisher_Zs_possible_vec = [(Fisher_Zs_min - 0.05*Fisher_Zs_gap):(Fisher_Zs_gap*0.01):(Fisher_Zs_max+0.05*Fisher_Zs_gap)];
            Fisher_Zs_possible_vec_size = size(Fisher_Zs_possible_vec)

            Fisher_Normal_fit_vec = (1/(sqrt(2*pi)*SAVE_FISHER_Z_STD{data_flag})) .* exp(-(Fisher_Zs_possible_vec-SAVE_FISHER_Z_MEAN{data_flag}).^2 ./ (2 .* SAVE_FISHER_Z_STD{data_flag}^2));
            Fisher_Normal_fit_vec_size = size(Fisher_Normal_fit_vec)

            SAVE_FISHER_NORMAL_FIT_VEC{data_flag} = Fisher_Normal_fit_vec;
            SAVE_FISHER_NORMAL_FIT_VEC_size = size(SAVE_FISHER_NORMAL_FIT_VEC)

        end

        if(rand_flag == mix_GAUSSIAN)
            num_of_iterations = 500;
            [Fisher_Zs_std,miu,prior, LogLike]=mixture_of_Gaussians(Fisher_Zs,num_of_Gaussians,num_of_iterations,Nsamples);



            % New ! Make a correction to the std. in order to
            % compensate for the noise of the samples
            % Note : There might be a problem if the STD is lower than what
            % we expect !!!
            Fisher_Zs_std
            std_correction = sqrt((1/(Nsamples-3)))
            Diff_Is = Fisher_Zs_std-sqrt((1/(Nsamples-3)))



            Fisher_Zs_size = size(Fisher_Zs)
            [hieght,bin_loc]=hist(Fisher_Zs,111);%bin_loc is Fisher_Zs_possible_vec
            Fisher_Zs_min = min(Fisher_Zs);  Fisher_Zs_max = max(Fisher_Zs); Fisher_Zs_gap = Fisher_Zs_max-Fisher_Zs_min;

            Fisher_Zs_possible_vec = [Fisher_Zs_min - 0.05*Fisher_Zs_gap:Fisher_Zs_gap*0.01:Fisher_Zs_max+0.05*Fisher_Zs_gap];
            Fisher_Zs_possible_vec_size = size(Fisher_Zs_possible_vec)
            clear y

            for i=1:num_of_Gaussians
                y(i,:)=prior(i)*1/(sqrt(2*pi)*Fisher_Zs_std(i))*exp(-(bin_loc-miu(i)).^2/(2*Fisher_Zs_std(i)^2));
            end

            Fisher_Normal_fit_vec = sum(y,1);%(bin_loc(2)-bin_loc(1))*5000
            Fisher_Normal_fit_vec_size = size(Fisher_Normal_fit_vec)
            SAVE_FISHER_NORMAL_FIT_VEC{data_flag} = Fisher_Normal_fit_vec;
            SAVE_FISHER_NORMAL_FIT_VEC_size = size(SAVE_FISHER_NORMAL_FIT_VEC)

            % Avoid fitting ! See what happens !
            Fisher_Zs_std = sqrt(Fisher_Zs_std.^2 -  (0/(Nsamples-3)))  % TRY FIT - TEMPORARY !!!!
            Fisher_Zs_std = sqrt(Fisher_Zs_std.^2 +  10*Z_sig_mean(end).^2./2);  % Here we take the sigma from the data and don't trust the 1/(n-3) !!!
        end

        %%     return;

        % Now calculate the desired N which is needed in order to get some
        % desired fraction
        i=1;
        %  return;

        % % % % %         [num_samples_needed(true_corr_flag+1) std_needed(true_corr_flag+1)] = ...
        % % % % % eps_vec = [ 0.5]; delta_vec = [ 0.2 ];
        % % % % % compute_num_samples_needed(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, eps_vec(1), delta_vec(1), Ngenes, miu, prior)
        % Now fill for a range of epsilon's and delta's
        %        return;
        %    eps_vec = [ 0.05 0.1 0.2 0.5 0.8]; delta_vec = [ 0.05 0.1 0.2 0.5 ];


        % % % % %         num_samples_table = zeros(length(eps_vec), length(delta_vec));
        % % % % %         for eps_ind = 1:length(eps_vec)
        % % % % %             for delta_ind = 1:length(delta_vec)
        % % % % %                 eps_ind
        % % % % %                 delta_ind
        % % % % %                 % Find n given eps,delta
        % % % % %                 [num_samples_table(eps_ind, delta_ind) std_dummy] = ...
        % % % % % compute_num_samples_needed(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, eps_vec(eps_ind), delta_vec(delta_ind), Ngenes, miu, prior);
        % % % % %
        % % % % %                 % Find delta given eps,n
        % % % % % %                delta_obtained = compute_confidence_delta(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, eps_vec(eps_ind), num_samples_table(eps_ind, delta_ind), Ngenes, miu, prior)
        % % % % % %                diff_deltas = delta_obtained-delta_vec(delta_ind)
        % % % % %                 % Find eps given delta,n
        % % % % % %                epsilon_obtained = compute_error_epsilon(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, delta_vec(delta_ind), num_samples_table(eps_ind, delta_ind), Ngenes, miu, prior)
        % % % % % %                diff_epsilons = epsilon_obtained-eps_vec(eps_ind)
        % % % % %             end
        % % % % %             eps_ind
        % % % % %         end


        % Now calculate different eps and delta for a given n, the number
        % of samples on the actual chip!!!
        % % % % %         eps_vec = [0.0001 0.001 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];  delta_vec = eps_vec;
        % % % % %         delta_obtained_vec = compute_confidence_delta(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, eps_vec, 78, Ngenes, miu, prior);
        % % % % %         eps_obtained_vec = compute_error_epsilon(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, delta_vec, 78, Ngenes, miu, prior);
        % % % % %         delta_obtained_vec
        % % % % %         eps_obtained_vec

        %            return;

        for(nsamples=nsamples_vec)

            done_nsamples = nsamples

            % Compute analytically the mean and (approx. ) variance
            % Here we cheat with the number of samples ....
           %%% dummy_nsamples = 3 + 1/Q_sig_best_fit;
            [inf_limit_frac(index, i) inf_limit_std(index, i) x_alpha ] = ...
                compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, 0.02 , nsamples-3, alpha,miu,prior);  % We use here the approximation of Fisher
                % In the above we're cheating. The '1' should be replaced
                % by Fisher_Zs_std
            
            % Now give the correct prob. of each gene to be in the top
            % N_top lists
            f_all_genes_inf_limit = zeros(2,Ngenes);
            f_all_genes_inf_limit(1,:) = [2 - normcdf((x_alpha*Fisher_Zs_std - Fisher_Zs).*sqrt(nsamples-3)) - normcdf((x_alpha*Fisher_Zs_std + Fisher_Zs).*sqrt(nsamples-3))].^2;  
            f_all_genes_inf_limit(2,:) = [normcdf((x_alpha*Fisher_Zs_std - Fisher_Zs).*sqrt(nsamples-3)) + normcdf((x_alpha*Fisher_Zs_std + Fisher_Zs).*sqrt(nsamples-3)) - 1].^2;   
            f_all_genes_inf_limit_ratio = 2 .* f_all_genes_inf_limit(1,:) ./ (1 + f_all_genes_inf_limit(1,:) + f_all_genes_inf_limit(2,:)); 
                
                
            
            if(data_flag  == RAND_DATA)

                % Now try a different randomization and see what changes
                f_all_genes_diff_limit = zeros(2,Ngenes);
                f_all_genes_diff_limit(1,:) = [2 - normcdf((x_alpha - Q_corr_vec).*sqrt(nsamples-3)) - normcdf((x_alpha + Q_corr_vec).*sqrt(nsamples-3))].^2;
                f_all_genes_diff_limit(2,:) = [normcdf((x_alpha - Q_corr_vec).*sqrt(nsamples-3)) + normcdf((x_alpha + Q_corr_vec).*sqrt(nsamples-3)) - 1].^2;
                f_all_genes_diff_limit_ratio = 2 .* f_all_genes_diff_limit(1,:) ./ (1 + f_all_genes_diff_limit(1,:) + f_all_genes_diff_limit(2,:));

                f_all_genes_diffdiff_limit = zeros(2,Ngenes);
                f_all_genes_diffdiff_limit(1,:) = [2 - normcdf((x_alpha - W_corr_vec).*sqrt(nsamples-3)) - normcdf((x_alpha + W_corr_vec).*sqrt(nsamples-3))].^2;
                f_all_genes_diffdiff_limit(2,:) = [normcdf((x_alpha - W_corr_vec).*sqrt(nsamples-3)) + normcdf((x_alpha + W_corr_vec).*sqrt(nsamples-3)) - 1].^2;
                f_all_genes_diffdiff_limit_ratio = 2 .* f_all_genes_diffdiff_limit(1,:) ./ (1 + f_all_genes_diffdiff_limit(1,:) + f_all_genes_diffdiff_limit(2,:));


            end

            
            %This works only for Gaussian, two-sided, true&sampled case:
            C_alpha = norminv(1-0.5*alpha);

            % Here we 'cheat' - calculate the correct asymptotic st.d. only
            % for the Gaussian case
            if(rand_flag == GAUSSIAN)
                asym_limit_std(index,i) = sqrt((1/sqrt(2*pi)) * exp(-C_alpha^2/2))/ (sqrt(2)*alpha*(pi*(nsamples-3)*Fisher_Zs_std^2)^(1/4)*sqrt(Ngenes));
            end

            % Now compute the std. 'stupidly', by assuming independence:
            stupid_independent_inf_limit_std(index, i) = sqrt(inf_limit_frac(index, i) * (1-inf_limit_frac(index, i)) / (alpha * Ngenes));

            kept_frac_ttt = cputime;

            % Now do also sampling (from model) and compare to analytic result
            [kept_frac_dist f_all_genes_from_data f_all_genes_from_data_one_NTOP] = sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, 0.02*Fisher_Zs_std, nsamples-3, Ngenes, alpha, iters, R, ...
                FALSE, miu, prior, 0.000000001 * Z_sig_std/Z_sig_mean, coeffs_vec, corr_mat);

            kept_frac_ttt = cputime - kept_frac_ttt
            samp_frac_mean(index, i) = mean(kept_frac_dist)/(alpha*Ngenes); samp_frac_std(index, i) = std(kept_frac_dist./(alpha*Ngenes));




            % Here we can approximate at least in the beginning the overlap from the data
            if((true_corr_flag == TWO_SAMPLED) & (nsamples<=Nsamples/2)) % For now skip this !!!!

                kept_frac_from_data_ttt = cputime;

                % Now do also sampling FROM THE DATA and compare to analytic result
                [ kept_frac_dist f_all_genes_from_data f_all_genes_from_data_one_NTOP]  = sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, 0.02*Fisher_Zs_std, nsamples, Ngenes, alpha, iters, R, ...
                    TRUE, miu , prior,  0.000000001 * Z_sig_std/Z_sig_mean, coeffs_vec, corr_mat);

                f_all_genes_from_data_vec = 2 .* f_all_genes_from_data(1,:) ./ (1 + f_all_genes_from_data(1,:) + f_all_genes_from_data(2,:)); 
                
                kept_frac_from_data_ttt = cputime - kept_frac_from_data_ttt
                samp_frac_mean_from_data(index, i) = mean(kept_frac_dist)/(alpha*Ngenes); samp_frac_std_from_data(index, i) = std(kept_frac_dist./(alpha*Ngenes));
            end

            % % % % % samp_frac_mean_from_data = samp_frac_mean; samp_frac_std_from_data = samp_frac_std; % Dummy wrong !!!!

            % Now do also sampling and compare to analytic result

            i=i+1;

        end  % loop on nsamples_vec


        % % % % %         % Now compute the whole distribution for the nsamples :
        % % % % %         % Now show as an example the bahaviour for one nsamples value
        % % % % %         % This is called only twice :
        % % % % %         if(true_corr_flag == TRUE_AND_SAMPLED) % i.e. true and sampled ..
        % % % % %             iters = 200; nsamples = 100; alpha  = alpha_vec(1); f_res = 0.025; soft_constrain_flag = 0;
        % % % % %             alpha
        % % % % %             CheckPermTopFunc(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, nsamples-3, Ngenes, alpha, iters, f_res, R, soft_constrain_flag, miu, prior);
        % % % % %         end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%  Compute once for 100, to comprare different datasets !!! 
        [inf_limit_frac_100 inf_limit_std_100 x_alpha_100 ] = ...
                compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, 100-3, alpha,miu,prior)  % We use here the approximation of Fisher
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        % adjust for Ngenes
        inf_limit_std(index,:) = inf_limit_std(index,:) ./ sqrt(Ngenes);

        if(rand_flag == mix_GAUSSIAN)
            asym_limit_std(index,:) = inf_limit_std(index, :);
        end

        % WRRRRONG !!!! JUST FOR NOW!!!!
        %%%%%%%%        samp_frac_std = inf_limit_std; samp_frac_mean = inf_limit_std;

        % Now plot the results :
        % Now print the figure. The simulation together with the delta of ngenes->infinity
        if(one_side_flag)
            side_str = 'one sided, ';
        else
            side_str = 'two sides, ';
        end
        if(true_corr_flag==TRUE_AND_SAMPLED)
            corr_str = 'true&samp., ';
        else
            corr_str = 'two samp., ';
        end
        if(rand_flag == GAUSSIAN)
            rand_str = 'gauss., ';
        end
        if(rand_flag == UNIFORM)
            rand_str = 'uni., ';
        end
        if(rand_flag == LINEAR)
            rand_str = 'lin., ';
        end
        if(rand_flag == mix_GAUSSIAN)
            rand_str = 'mix. gauss., ';
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure; hold on; % Make new figure not to distroy the old one !!!!
        % subplot(2,2,true_corr_flag+1); hold on;

        % Analytic saddle-point approximation , here 'easy' stds.
        plot(nsamples_vec, inf_limit_frac(index,:) - stupid_independent_inf_limit_std(index,:), 'g+');
        plot(nsamples_vec, inf_limit_frac(index,:) + stupid_independent_inf_limit_std(index,:), 'g+');

        % Sampled from the Gaussian model
        errorbar(nsamples_vec, samp_frac_mean(index,:), samp_frac_std(index,:), 'r');

        if(true_corr_flag == TRUE_AND_SAMPLED)
            % Analytic saddle-point approximation , here correct stds.
            errorbar(nsamples_vec, inf_limit_frac(index,:), inf_limit_std(index,:));

            % Just plot the alpha, which is the worst possible case
            plot(nsamples_vec, repmat(alpha, 1, length(nsamples_vec)), ':c');

            %%%legend('ind-saddle_{minus}', 'ind-saddle_{plus}',
            legend('ind. std.', 'ind. std.',  'sampled from model', 'Saddle approx.', '\alpha (RAND)', 2);
            % % %             legend('sampled from model','Saddle approx.', 2);
        else
            % Give it up for now , bad fit :
            %            errorbar(nsamples_vec(1:small_vec_len), samp_frac_mean_from_data(index,1:small_vec_len), samp_frac_std_from_data(index,1:small_vec_len), 'm');

            plot(nsamples_vec(1:small_vec_len), samp_frac_mean_from_data(index,1:small_vec_len), 'm*');

            % Analytic saddle-point approximation , here correct stds.
            errorbar(nsamples_vec, inf_limit_frac(index,:), inf_limit_std(index,:));
            %%%%            plot(nsamples_vec, inf_limit_frac(index,:));
            %%%   legend('saddle_{minus}', 'saddle_{plus}',

            % Just plot the alpha, which is the worst possible case
            plot(nsamples_vec, repmat(alpha, 1, length(nsamples_vec)), ':c');

            legend( 'ind. std.', 'ind. std.',  'sampled from model','sampled from data','Saddle approx.', '\alpha (RAND)',  2);  % just before last, add :  'sampled from data',
            % % %             legend('sampled from model', 'Saddle approx.', 2);
        end


        
        
        ylabel('Mean kept Frac f*'); xlabel('No. Samples. ');
        %        xlabel(['No. Samples. The number needed for ' num2str(1-epsilon_needed) ' with confidence ' num2str(1-delta_needed) ' is: '  num2str(num_samples_needed(true_corr_flag+1))]);
        
        SS_Err = sum((inf_limit_frac(index,:) - samp_frac_mean_from_data(index,1:small_vec_len)).^2 ./ inf_limit_frac(index,:)) 
        title([data_str ' Data, f dist. ' rand_str corr_str side_str 'Ngenes=' num2str(Ngenes) ' alpha=' num2str(alpha) ' N_{TOP}= ' num2str(round(Ngenes*alpha)) ' iters= ' num2str(iters)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        % Now plot only the std. !!!
        figure; hold on; title([ 'Only STDs !!!'  corr_str]);
        plot(nsamples_vec, stupid_independent_inf_limit_std(index,:), 'g+');
        plot(nsamples_vec, samp_frac_std(index,:), 'r');
        plot(nsamples_vec, inf_limit_std(index,:));
        plot(nsamples_vec, asym_limit_std(index,:), 'm');
        legend('ind. std.', 'sampled from model', 'Saddle approx.', 'Asym approx.', 2);


        % % % %         % Now plot only the std. Ratio !!!
        % % % %         figure; hold on; title([ 'STDs Ratio !!!'  corr_str]);
        % % % %         plot(nsamples_vec, asym_limit_std(index,:)./inf_limit_std(index,:), 'm');
        % % % %         legend('Ratio', 2);

        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        index=1;
        figure; hold on; subplot(2,2,1); hold on; plot(f_all_genes_inf_limit(1,:), f_all_genes_from_data(1,:), '.'); 
        title(['Comparison of analytic and data simulation for each gene, n_s = ' num2str(max(nsamples_vec)) ' data is ' data_str]);
        xlabel('Analytic'); ylabel('Data'); plot([0:0.01:1],[0:0.01:1], 'r'); 
        
        subplot(2,2,2); hold on; plot(Fisher_Zs, f_all_genes_from_data(1,:), '.');   plot(Fisher_Zs, f_all_genes_inf_limit(1,:), '.r'); 
        plot(Fisher_Zs, repmat(mean(f_all_genes_from_data(1,:)./alpha), 1, Ngenes), 'g'); 
        plot(Fisher_Zs, repmat(mean(f_all_genes_inf_limit(1,:)./alpha), 1, Ngenes), 'm'); 
        plot(Fisher_Zs, repmat(samp_frac_mean_from_data(index,end), 1, Ngenes), 'x');
        plot(Fisher_Zs, repmat(inf_limit_frac(index,end), 1, Ngenes), 'xr');

        plot(Fisher_Zs, f_all_genes_from_data_one_NTOP(1,end), 'c.'); 
        plot(Fisher_Zs, repmat(mean(f_all_genes_from_data_one_NTOP(1,end)./alpha), 1, Ngenes), 'c+'); 
        
        if(data_flag == RAND_DATA)
            plot(Q_corr_vec, f_all_genes_diff_limit(1,end), 'k.');
            plot(Q_corr_vec, repmat(mean(f_all_genes_diff_limit(1,end)./alpha), 1, Ngenes), 'k+');
        end
        
        title('Overlap for each gene'); xlabel('Fisher Z'); ylabel('Prob. in N_{TOP}'); 
        legend('From Data', 'Analytic', 'Mean Data', 'Mean Analytic', 'Mean f data', 'Mean f Analytic', 'From Data Sqr', 'Mean From Data Sqr', 'Analytic diff', 'Mean Analytic diff');
        
        
        
         subplot(2,2,3); hold on;         
        
         if(data_flag == RAND_DATA)
             plot(Q_corr_vec, f_all_genes_diff_limit(1,end), 'k.');
             plot(Q_corr_vec, repmat(mean(f_all_genes_diff_limit(1,end)./alpha), 1, Ngenes), 'k+');

             plot(W_corr_vec, f_all_genes_diffdiff_limit(1,end), 'm.');
             plot(W_corr_vec, repmat(mean(f_all_genes_diffdiff_limit(1,end)./alpha), 1, Ngenes), 'm+');
         end
         plot(Fisher_Zs, repmat(inf_limit_frac(index,end), 1, Ngenes), 'r');
        
        title('Overlap for each gene'); xlabel('Fisher Z'); ylabel('Prob. in N_{TOP}'); 
        legend('From Data', 'Analytic', 'Mean Data', 'Mean Analytic', 'Mean f Analytic');
        
%         subplot(2,2,3); hold on; plot(f_all_genes_inf_limit(1,:), f_all_genes_from_data_one_NTOP(1,:), '.'); 
%         title(['Comparison of analytic and data SQR simulation for each gene, n_s = ' num2str(max(nsamples_vec)) ' data is ' data_str]);
%         xlabel('Analytic'); ylabel('Data SQR'); plot([0:0.01:1],[0:0.01:1], 'r'); 

        
        
        
        subplot(2,2,4); hold on; plot(f_all_genes_from_data(1,:), f_all_genes_from_data_one_NTOP(1,:), '.'); 
        title(['Comparison of data and data SQR simulation for each gene, n_s = ' num2str(max(nsamples_vec)) ' data is ' data_str]);
        xlabel('Data'); ylabel('Data SQR'); plot([0:0.01:1],[0:0.01:1], 'r');                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 5.5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% % % 
% % %         % Now plot many Gaussians, according to Eytan suggestion
% % %         figure; hold on;
% % %         color_vec = 'brgmkcy:o*+sdxv<>^';
% % %         f_res = 0.001;
% % %         f_vec = [f_res:f_res:1-f_res];
% % %         for color_ind = 1:length(nsamples_vec)
% % %             %            plot(f_vec, (sqrt(Ngenes)/(sqrt(2*pi) * inf_limit_std(index,color_ind)  )), color_vec(color_ind));
% % %             plot(f_vec, (1/(sqrt(2*pi) * inf_limit_std(index,color_ind)  )) .* ...
% % %                 exp ( -(f_vec -  inf_limit_frac(index,color_ind)).^2 .* 1 ./(2.0 * inf_limit_std(index,color_ind).^2)) , color_vec(color_ind));
% % %         end
% % % 
% % %         legend([num2str(nsamples_vec(1)) ' samples'], [num2str(nsamples_vec(2)) ' samples'], [num2str(nsamples_vec(3)) ' samples'], ...
% % %             [num2str(nsamples_vec(4)) ' samples'], [num2str(nsamples_vec(5)) ' samples']);
% % %         title('Prob. density func. of f for different number of samples'); xlabel('f'); ylabel('Pr(f)');
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
        
        
        

        % %      Do sub-sampling. Currently we try to avoid it.
        % %         else  % Here do some sampling-based method
        % %
        % %
        % %             corr_num_iters = 10000; % number of iterations to randomize couples
        % %
        % %
        % %             for iter = 1:corr_num_iters
        % %                 % randomize the couples
        % %                 rand_perm = randperm(1,Nsamples);
        % %
        % %                 for i=1:2:floor(Nsamples/2)*2-1
        % %                     % Calculate the current two-point data and normalize it
        % %                     cur_data_vec = R.data([rand_perm(i), rand_perm(i+1)],:);
        % %                     cur_data_vec = cur_data_vec - repmat(mean(cur_data_vec, 2), 1, 2);
        % %                     cur_data_vec = cur_data_vec ./ repmat( sqrt(sum( cur_data_vec .^2, 1)), 1, 2);
        % %
        % %                     cur_labels_vec = R.labels(rand_perm(i), rand_perm(i+1));
        % %
        % %                     % Calculate the correlations
        % %                     curr_corr_vec = cur_data_vec .* repmat(cur_labels_vec, 1, Ngenes);
        % %
        % %                     for j=1:Ngenes % for now do correlation gene by gene !
        % %
        % %                         garbage  = curr_corr_vec;
        % %
        % %                     end
        % %
        % %                 end
        % %             end
        % %         end % if Fisher flag

        plot_ind = plot_ind+1;

        index = index + 1;

        %%%      end  % loop on alpha
        %%% end  % loop on two samps flag


    end  % loop on true corrs

end % loop on iters_vec
    
    % Now plot the Fisher Z's for comparison
    % % % % % %     SAVE_FISHER_Z{data_flag} = min(SAVE_FISHER_Z{data_flag}, 0.31);
    % % % % % %     SAVE_FISHER_Z{data_flag} = max(SAVE_FISHER_Z{data_flag}, -0.345);
    % % % % % %     SAVE_FISHER_NORMAL_FIT_VEC{data_flag} = ...
    % % % % % %         (1/(sqrt(2*pi)*SAVE_FISHER_STD{data_flag})) .* exp(-(Fisher_Zs_possible_vec-SAVE_FISHER_MEAN{data_flag}).^2 ./ (2 .* SAVE_FISHER_STD{data_flag}^2));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the q_z distribution
    figure; %subplot(2,2,3);
    hold on;  num_bins = max(250, floor(Ngenes/250)); % % num. of bins is arbitrary. Should be determined more cleverly
    [H bin_pos] = hist(SAVE_FISHER_Z{data_flag}, num_bins);
    hist(SAVE_FISHER_Z{data_flag}, num_bins); xlabel('Correlation'); ylabel('Freq.'); title(['Data is: ' data_str ' Hist. of q_z, all genes Fisher Zs']);
    plot(Fisher_Zs_possible_vec, SAVE_FISHER_NORMAL_FIT_VEC{data_flag} .* Ngenes .* (bin_pos(2)-bin_pos(1)), 'r');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set now wider axis
    AX = gca;
    X_ax_lim = get(AX, 'xlim'); Y_ax_lim = get(AX, 'ylim');
    X_ax_lim(1) = min(X_ax_lim(1), X_ax_lim(1)*1.2);     X_ax_lim(2) = max(X_ax_lim(2), X_ax_lim(2)*1.2);
    Y_ax_lim(1) = min(Y_ax_lim(1), Y_ax_lim(1)*1.2);     Y_ax_lim(2) = max(Y_ax_lim(2), Y_ax_lim(2)*1.2);
    set(AX, 'xlim', X_ax_lim); set(AX, 'ylim', Y_ax_lim);



    if(rand_flag == GAUSSIAN)
        data_mean = mean(SAVE_FISHER_Z{data_flag});
        data_std = std(SAVE_FISHER_Z{data_flag});
        saved_mean = SAVE_FISHER_Z_MEAN{data_flag};
        saved_std = SAVE_FISHER_Z_STD{data_flag};

        % Now take a random sample from the same Gaussian distribution and see
        % if it looks the same
        rand_gauss_samp = SAVE_FISHER_Z_STD{data_flag} .* randn(1, Ngenes) + SAVE_FISHER_Z_MEAN{data_flag};
        figure; % subplot(2,2,4);
        hold on;  num_bins = max(250, floor(Ngenes/250)); % % num. of bins is arbitrary. Should be determined more cleverly
        [H bin_pos] = hist(rand_gauss_samp, num_bins);
        hist(rand_gauss_samp, num_bins); xlabel('Correlation'); ylabel('Freq.'); title([data_str ' Data Hist. of rand. samp.']);
        plot(Fisher_Zs_possible_vec, SAVE_FISHER_NORMAL_FIT_VEC{data_flag} .* Ngenes .* (bin_pos(2)-bin_pos(1)), 'r');
        AX = gca;
        set(AX, 'xlim', X_ax_lim); set(AX, 'ylim', Y_ax_lim);
    end


end % end loop on data types
% cd data
% save OR_WANG_RESULTS nsamples_vec inf_limit_frac inf_limit_std samp_frac_mean samp_frac_std samp_frac_mean_from_data samp_frac_std_from_data





%%%%%%%%%%%%%%%%%%%%% Example of doing Gaussian fit
% x_vec = [-10:0.05:10]; Z_MEAN = 3; Z_STD = 3;
% Normal_fit_vec = (1/(sqrt(2*pi)*Z_STD)) .* exp(-(x_vec-Z_MEAN).^2 ./ (2 .* Z_STD^2));
% Ngenes = 10000; rand_gauss_samp = Z_STD .* randn(1, Ngenes) + Z_MEAN;
% num_bins = max(250, floor(Ngenes/250)); % % num. of bins is arbitrary. Should be determined cleverly
% [H bin_pos] = hist(rand_gauss_samp, num_bins);
% figure; hold on; hist(rand_gauss_samp, num_bins); xlabel('Value'); ylabel('Freq.'); title([' Data Hist. of rand. samp.']);
% plot(x_vec, Normal_fit_vec.* Ngenes .* (bin_pos(2)-bin_pos(1)), 'r');
%%%%%%%%%%%%%%%%%%%%%


time_elapsed = cputime -ttt_time










