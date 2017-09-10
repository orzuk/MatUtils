% Here we load a data, compute the correlations and get the desired
% fraction
%  path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
%  path(path, 'E:\Research\PACPerm\numeric');
% path(path, 'C:\Weizmann\Research\PACPerm\numeric');
% path(path, 'C:\Weizmann\Research\PACPerm\numeric\nmm\linalg');
% 
% path(path, '/a/fs-01/pdoc/mrosenz/Liat/GeneCorr/Or');
path(path, 'E:\Gene_Corr\');
TRUE = 1; FALSE = 0;


kept_fraction=0.01;
KNN = 10; % parameter for knn missing values algorithm
Fisher_R_to_Z = 1; % Flag saying if to perform fisher RtoZ transformation
remove_zero_flag = FALSE;   % Flag saying to remove all the features with correlation
% zero from the data, since they cause bias.


% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; FROM_DATA = 3;mix_GAUSSIAN=4;
student_t=5;% in the last one we take q simply according to the data bins
rand_flag = GAUSSIAN;

% data files
OLD_VANT_VEER = 1; NEW_ROSSETA = 2; WANG = 3;     % gene expression
TOPIC = 4; % Author-topic matching
NIPS_DOROTHEA = 5; NIPS_ARCENE = 6; NIPS_DEXTER = 7; digits = 8; % Various datasetes from NIPS contest
RAND_DATA = 9; ROSENWALD=10;HC=11;% Here we randomize a data so that we have control on it, and also can work with matlab 6.5 (no loading from files)
Bhattacharjee=12;BEER=13;YEOH=14;

% Vector saying which datasets come in a sparse form
IS_SPARSE_VEC = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];
CALC_CORRS_VEC = [1,1,1,0,0,0,0,1,1,1,1,1,1,1];  % calc for Gizette [1,1,1,0,0,0,0,0];
CALC_VAR_VEC = [1,0,1,0,0,0,0,0,0,0,0,0,0,0];

calc_corrs_flag = FALSE; % flag saying if to calculate the correlations from the data, or do something else

% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TWO_SAMPLED;


%alpha = 0.014; 
%TOP_GENES = 70;
a=exp(1); b=-1;
c=a; d=b;
alpha=[];
miu=[];
prior=1;
%nsamples=3;
iters = 20;
alpha_vec =  [ 0.012  0.12]; % corresponds to ~70 and ~700 genes





% Now chose which data to load - New : Run on all data's !!!!


for data_flag =  1:3 OLD_VANT_VEER % WANG
 % NEW_ROSSETA % digits; %NIPS_DEXTER; %%NIPS_ARCENE; %NIPS_DOROTHEA; %
 % TOPIC; % WANG; %%% NEW_ROSSETA; %%%OLD_VANT_VEER;%%ROSENWALD
    cur_data_is = data_flag%%%%%%HC%%%Bhattacharjee%%%%BEER%%%YEOH
    Fisher_R_to_Z = 1; % Flag saying if to perform fisher RtoZ transformation
    one_side_flag = TWO_SIDES; % make sure we do two-sided unless later we change the flag
    calc_corrs_flag = CALC_CORRS_VEC(data_flag);
    estimate_variance_flag = CALC_VAR_VEC(data_flag);
    % Now for each data we do one figure !
    plot_ind=1;
    figure; subplot(2,2,1); hold on;


    clear R; % free memory
    if(data_flag == HC)
        rand_flag = GAUSSIAN;
%         num_of_Gaussians=2;
%         num_of_itterations=100;
        kept_fraction=0.5;
        d_var_n=0;
         res = 4; % the resulution of N that we use
        max_nsamples = 30; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
       load('../data/HC'); % The name of the data file containing the expression and labels.
        data_str = 'HC';
        data(data<1)=1;
        data=log2(data);
        present=real(pval<0.05);
        present_idx=find(sum(present,2)>=size(pval,2)*kept_fraction);
        present_data=data(present_idx,:);
        R.dat=present_data;
        R.Labels=labels';
          b = -1.1999;
  a = exp(0.76951);

        b = -1.2093;
        a = exp(0.79232);
        c = 3.885;
        d = -1.442;
%         a=c;
%         b=d;
%  b = -1.7466
%   a = 1.8441
    end
    if(data_flag == BEER)
        rand_flag = GAUSSIAN;
%         num_of_Gaussians=2;
%         num_of_itterations=100;
        kept_fraction=0.5;
        d_var_n=0;
         res = 4; % the resulution of N that we use
        max_nsamples = 43; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
       load('../data/beer'); % The name of the data file containing the expression and labels.
        data_str = 'BEER';
        pval=data<32;
        data(data<32)=32;
        data=log2(data);
        present=real(pval<0.05);
        present_idx=find(sum(present,2)>=size(pval,2)*kept_fraction);
        present_data=data(present_idx,:);
        R.dat=present_data;
        R.Labels=labels;
        b = -1.1594;
        a = exp(0.64195);
        c = 2.873; 
        d = -1.325;
a=c
b=d
    end

    if(data_flag == YEOH)
        rand_flag = GAUSSIAN;
%         num_of_Gaussians=2;
%         num_of_itterations=100;
        kept_fraction=.5;
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
%         Y = var(data,0,2);
%         [val,ind]=sort(-Y);
%         present_data=data(ind(1:4000),:);
        present_idx=find(sum(present,2)>=N_s*kept_fraction);
        present_data=data(present_idx,:);
        R.dat=present_data;
        R.Labels=labels;
%       p1 = -1.027;
%       p2 = 7.0523;
%   b = -1.0923
%  a = 7.5463
% %  b = -1.7466
%   a = 1.8441
    end


    if(data_flag == Bhattacharjee)
        rand_flag = GAUSSIAN;
%         rand_flag = mix_GAUSSIAN;
%         num_of_Gaussians=2;
%         num_of_itterations=100;
        kept_fraction=0.5;
        d_var_n=0;
%         remove_zero_flag = TRUE;
         res = 4; % the resulution of N that we use
        max_nsamples = 33; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        min_nsamples = 9; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
       load('../data/Bhattacharjee_66.mat'); % The name of the data file containing the expression and labels.
        data_str = 'Bhattacharjee';
        data(data<30)=30;
        data=log2(data);
        Y = var(data,0,2);
        [val,ind]=sort(-Y);
        present_data=data(ind(1:4000),:);
        R.dat=present_data;
        R.Labels=labels;
        b = -1.2408;
        a = exp(0.89492);

       c =       8.583 ;
       d =       -1.84;
       a=c
       b=d
%        c =     0.03471 var of true dist
%           b = -1.6195;
%           a = 1.5354;

    end


    if(data_flag == OLD_VANT_VEER)
        rand_flag = GAUSSIAN;
         res = 4; % the resulution of N that we use
        max_nsamples = 48; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
       load('../data/breastData.mat'); % The name of the data file containing the expression and labels.
        data_str = 'VantVeerOld';
        d_var_n=0;
      b = -1.1847;
      a = exp(0.7479);
      c =2.698 ;
      d =-1.284;
      a=c
      b=d
    end
    if(data_flag == NEW_ROSSETA)
         res = 10; % the resulution of N that we use
        max_nsamples = 100; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
        loading_big_rosseta = 1;
        load('../data/Rosetta_data.mat'); % The name of the data file containing the expression and labels.
        data_str = 'VantVeerNew';
        rand_flag = GAUSSIAN;
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
        
        R.Labels = real(Labels<100); R.Labels = R.Labels';
        loading_big_rosseta = 1
        d_var_n=-3;
    end

    if(data_flag == WANG)
        rand_flag = GAUSSIAN;
        load('../data/WANG_DATA.mat');
        data_str = 'Wang';
        kept_fraction=0;
         res = 10; % the resulution of N that we use
        max_nsamples = 100; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
        gene_expression_log2(gene_expression_log2<5)=5;
        Y = var(gene_expression_log2,0,2);
        [val,ind]=sort(-Y);
        present_data=gene_expression_log2(ind(1:4000),:);
        R.dat=present_data;
        R.dat = R.dat + 0.001*rand(size(R.dat,1),size(R.dat,2));
%%%%%%%removing genes with too many missing values
%         present=real(pval<0.05);
%         present_idx=find(sum(present,2)>=size(pval,2)*kept_fraction);
%         present_data=gene_expression_log2(present_idx,:);
%         data=present_data;
%         % Here we must do log since Gary didn't do it
%        R.dat = data;;
        R.Labels = Prognosis';
        calc_corrs_flag = TRUE;
        d_var_n=-3;
    end

    if(data_flag == ROSENWALD)
        rand_flag = GAUSSIAN;
%         num_of_Gaussians=1;
%         num_of_itterations=100;
%        a=1;
%         b=-1;
        load('../data/ROSEN1.mat');
        data_str = 'ROSENWALD';
         kept_fraction=.5
         res = 10; % the resulution of N that we use
        max_nsamples = 120; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
%%%%%%%removing genes with too many missing values
%         Here we must do log since Gary didn't do it
        present=real(pval<0.05);
        present_idx=find(sum(present,2)>=size(pval,2)*kept_fraction);
        present_data=data(present_idx,:);
        R.dat=present_data;
        R.Labels = labels;
        calc_corrs_flag = TRUE;
        d_var_n=0;
        b = -1.0961
        a = exp(0.43978);
        c =       2.476; 
        d =      -1.283 ;
        a=c;
        b=d;
%        c =    0.009813  (0.009461, 0.01017)


    end
    if(data_flag == TOPIC)
        load('../data/topic_data.mat');
        data_str = 'Author-Topic';
        selected_idx=find(sum(Topics_N_counts,2)~=0);
        R.dat = Topics_N_counts(selected_idx,:);
        R.Labels = N_Labels';
        calc_corrs_flag = TRUE;   % Many zero correlations in this dataset.
        rand_flag = mix_GAUSSIAN;
        num_of_Gaussians=2;
        num_of_itterations=100;
        d_var_n=0;
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

    if(data_flag == digits)
        load('../data/digits.mat');
        data_str = 'NIPS-Gisette';
        rand_flag = mix_GAUSSIAN;
        num_of_Gaussians=3;
        num_of_itterations=100;
%         rand_flag =student_t;
        R.dat = data;
        
        % avoid constant rows by randomizing them
        const_rows_indexes = find(max(R.dat,[],2) == min(R.dat,[],2));
        R.dat(const_rows_indexes,:) = rand(length(const_rows_indexes), size(R.dat,2));
        res = 10; % the resulution of N that we use
        min_nsamples = 10; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        max_nsamples = 200; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)
        nsamples_vec = min_nsamples:res:max_nsamples;
        
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
        
        R.Labels = samples';
        calc_corrs_flag = TRUE;
%         [val,indd]=sort(var_corr(:,10));
%         R.dat=R.dat(indd(300:600),:);

%         Pearson_corrs=Pearson_corrs(indd(2000:end));
        d_var_n=0;
    end


    if(data_flag == RAND_DATA)
        
        Ngenes = 5000; Nsamples = 100; 
        R.dat = rand(Ngenes,Nsamples);
        R.Labels = round(rand(1, Nsamples)); % 0/1 labels 
        
    end
    

    Ngenes = size(R.dat, 1)
    Nsamples = size(R.dat, 2)



    if(remove_zero_flag == TRUE)  % remove zero correlations only if we do not calculate them
        non_zero_ind = find(sum(R.dat==0,2)~=size(R.dat,2));
        R.dat = R.dat(non_zero_ind,:);
        Ngenes = size(R.dat, 1);
        if(calc_corrs_flag ~= TRUE)
            Pearson_corrs=Pearson_corrs(non_zero_ind);
        end

    end

    if(calc_corrs_flag == TRUE)
        % Calculate the mean correlation of each gene with survival
        normed_labels = R.Labels - mean(R.Labels); normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
        normed_data = R.dat;
%         normed_data(4965,:) = rand(1,Nsamples);
        normed_data = normed_data - repmat(mean(normed_data, 2), 1, Nsamples);
        normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, Nsamples);
        %normed_data = normed_data';

        Pearson_corrs = sum(normed_data .* repmat(normed_labels, Ngenes, 1),2);
        size(Pearson_corrs)
    % here the correlations are already given
    end

    Ngenes = size(R.dat, 1)

    inf_limit_frac = zeros(2 * length(alpha_vec), length(nsamples_vec));
    inf_limit_std = zeros(2 * length(alpha_vec), length(nsamples_vec));
    samp_frac_mean = zeros(2 * length(alpha_vec), length(nsamples_vec));
    samp_frac_std = zeros(2 * length(alpha_vec), length(nsamples_vec));
    samp_frac_mean_from_data = zeros(2 * length(alpha_vec), length(nsamples_vec));
    samp_frac_std_from_data = zeros(2 * length(alpha_vec), length(nsamples_vec));

        % Now try to estimate the variance of the correlation of each gene with
        % survival. The method we use is go over all the couples.
        % Below the 2-factor is for both values of true_corr_flag


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

        SAVE_FISHER_Z{data_flag} = Fisher_Zs; % Save the Fisher Z's for plotting later

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
            
            
            % New ! Make a correction to the std. in order to
            % compensate for the noise of the samples
            Fisher_Zs_std = sqrt(Fisher_Zs_std^2 -  c*(Nsamples+d_var_n)^d);  % TRY FIT - TEMPORARY !!!!

            Fisher_Zs_min = min(Fisher_Zs);  Fisher_Zs_max = max(Fisher_Zs); Fisher_Zs_gap = Fisher_Zs_max-Fisher_Zs_min;

            Fisher_Zs_possible_vec = [Fisher_Zs_min - 0.05*Fisher_Zs_gap:Fisher_Zs_gap*0.01:Fisher_Zs_max+0.05*Fisher_Zs_gap];
            Fisher_Normal_fit_vec = (1/(sqrt(2*pi)*SAVE_FISHER_Z_STD{data_flag})) .* exp(-(Fisher_Zs_possible_vec-SAVE_FISHER_Z_MEAN{data_flag}).^2 ./ (2 .* SAVE_FISHER_Z_STD{data_flag}^2));
            SAVE_FISHER_NORMAL_FIT_VEC{data_flag} = Fisher_Normal_fit_vec;
            

        end

        if(rand_flag == mix_GAUSSIAN)
            % New ! Make a correction to the std. in order to
            % compensate for the noise of the samples
            [Fisher_Zs_std,miu,prior]=mixture_of_Gaussians(Fisher_Zs,num_of_Gaussians,num_of_itterations,Nsamples);
            Fisher_Zs_std = sqrt(Fisher_Zs_std.^2 -  (c*(Nsamples+d_var_n)^d));  % TRY FIT - TEMPORARY !!!!

            [hieght,bin_loc]=hist(Fisher_Zs,100);%bin_loc is Fisher_Zs_possible_vec
            Fisher_Zs_min = min(Fisher_Zs);  Fisher_Zs_max = max(Fisher_Zs); Fisher_Zs_gap = Fisher_Zs_max-Fisher_Zs_min;

            Fisher_Zs_possible_vec = [Fisher_Zs_min - 0.05*Fisher_Zs_gap:Fisher_Zs_gap*0.01:Fisher_Zs_max+0.05*Fisher_Zs_gap];

            clear y
           
            for i=1:num_of_Gaussians
                y(i,:)=prior(i)*1/(sqrt(2*pi)*Fisher_Zs_std(i))*exp(-(bin_loc-miu(i)).^2/(2*Fisher_Zs_std(i)^2));
            end

            Fisher_Normal_fit_vec = sum(y,1);%(bin_loc(2)-bin_loc(1))*5000
            SAVE_FISHER_NORMAL_FIT_VEC{data_flag} = Fisher_Normal_fit_vec;
            

        end
%          epsilon_needed=0.7; delta_needed=0.1;true_corr_flag=1;
%        [num_samples_needed(true_corr_flag+1) std_needed(true_corr_flag+1)] = compute_num_samples_needed(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, epsilon_needed, delta_needed, Ngenes)
%         delta_obtained = compute_confidence_delta(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, epsilon_needed, 78, Ngenes)
%         epsilon_obtained = compute_error_epsilon(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, alpha, delta_needed, 78, Ngenes)

        index = 1;
        for true_corr_flag=0:0

        % Now calculate the desired N which is needed in order to get some
        % desired fraction
      rescaled_fisher_ZS=(Fisher_Zs-Fisher_Zs_mean)*Fisher_Zs_std/sqrt(Fisher_Zs_std.^2 + (c*(Nsamples+d_var_n)^d))+Fisher_Zs_mean;
       i=1;
   
       
       
       
       
       if( (calc_corrs_flag == TRUE) & (estimate_variance_flag == TRUE))
           corrs_iters = 100;
           [N_VEC Pearson_corrs Fisher_Zs Rho_sig_mean Rho_sig_std Rho_bias_mean Rho_bias_std Z_sig_mean Z_sig_std Z_bias_mean ...
               Z_bias_std chunk_Rho_sig chunk_Z_sig chunk_Rho_bias chunk_Z_bias ] = ...
               CalcZVarianceFromDataFunc(R, alpha_vec, rand_flag, corrs_iters ); % take only max samples to save time
           subplot(2,3,4); hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. N=' num2str(N_VEC)]);
           xlabel('true P'); ylabel('std');
           subplot(2,3,5); hold on; plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. N='  num2str(N_VEC(1))]);
           xlabel('true Z'); ylabel('std');



           corrs_iters2=100;


           % Skip the following: Very time-consuming
           [N_VEC2 N_HALF_VEC2 Z_sig_mean_plot Q_sig_mean_plot ] = CalcTwoZVariancesFromDataFunc(R, alpha, rand_flag, corrs_iters2 );

           

             

           % Now do fit
           linfit = fittype('a*x+b');
           Z_fit = fit(log(N_HALF_VEC2'), log(Z_sig_mean_plot'), linfit);
           a = Z_fit.a; b = Z_fit.b;
           % Note: This gives different (and crappy) results ! 
           orig_fit = fittype('a*x^b');
           ZZ_fit = fit(N_HALF_VEC2', Z_sig_mean_plot', orig_fit);
           comp_fit = fittype('f+c*x^(-0.4)');
           Q_fit = fit(N_VEC2', Q_sig_mean_plot', comp_fit);
           
           
           % Plot and compare both variances
           figure; hold on; plot(N_HALF_VEC2, Z_sig_mean_plot, '.'); plot(N_VEC2, Q_sig_mean_plot, '.r');  plot(N_HALF_VEC2, Q_sig_mean_plot(1:length(N_HALF_VEC2))-Z_sig_mean_plot, '.g');
           plot(N_HALF_VEC2, 1 ./ sqrt(N_HALF_VEC2-3), 'k');    plot(N_VEC2, std(Fisher_Zs) + 1 ./ sqrt(N_VEC2-3), 'k');
           plot(N_HALF_VEC2,  log(N_HALF_VEC2).^Z_fit.a .* Z_fit.b, 'm');
           legend('Ave Z std.', 'Q std.', 'diff', 'Fisher', 'Fisher plus one'); xlabel('Nsamples'); ylabel('Std');
           title(['Std. of each genes Z and of total Q, Data is ' data_str]);

           
            % Plot and compare both variances on LogLog scale
           figure; hold on; plot(log(N_HALF_VEC2), log(Z_sig_mean_plot), '.'); 
           plot(log(N_VEC2), log(Q_sig_mean_plot), '.r');  
           plot(log(N_HALF_VEC2), Z_fit.a .* log(N_HALF_VEC2) + Z_fit.b, 'm');
           
           legend('Log Ave Z std.', 'Log Q std.', 'Log Z fit');  xlabel('Log Nsamples'); ylabel('Log Std');
           title(['LogLog plot of Std. of each genes Z and of total Q, Data is ' data_str]);
         
            
                    
           
       end   % plot of std of the

       
       
        return;
       
       
       for(nsamples=nsamples_vec)
            i

            done_nsamples = nsamples
            nsam_eff=1/(a*(nsamples+d_var_n)^b);
            nsamples_dist=1/(c*(nsamples+d_var_n)^d);
            % Compute analytically the mean and (approx. ) variance
            [inf_limit_frac(index, i) inf_limit_std(index, i) x_alpha ] = compute_inf_limit_fraction_type2(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std,nsam_eff,nsamples_dist ,alpha,miu,prior);  % We use here the approximation of Fisher
%             [inf_limit_frac(index, i) inf_limit_std(index, i) x_alpha ] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std,nsam_eff ,alpha,miu,prior);  % We use here the approximation of Fisher
            if(i==length(nsamples_vec))
                saved_x_alpha=x_alpha;
            end
            % Now compute the std. 'stupidly', by assuming independence:
            stupid_independent_inf_limit_std(index, i) = sqrt(inf_limit_frac(index, i) * (1-inf_limit_frac(index, i)) / (alpha * Ngenes));
            % Now do also sampling and compare to analytic result
%             kept_frac_dist = sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, 1/(a*(nsamples+d_var_n)^b), Ngenes, alpha, iters, R, FALSE,miu,prior);
%             samp_frac_mean(index, i) = mean(kept_frac_dist)/(alpha*Ngenes); samp_frac_std(index, i) = std(kept_frac_dist./(alpha*Ngenes));
            data_kept_frac_dist=sample_kept_fraction_data_distribution(rand_flag, one_side_flag, true_corr_flag, rescaled_fisher_ZS', nsam_eff, Ngenes, alpha, iters, R, FALSE);
            data_samp_frac_mean(index, i) = mean(data_kept_frac_dist)/(alpha*Ngenes);data_samp_frac_std(index, i) = std(data_kept_frac_dist./(alpha*Ngenes));
    
            % Here we can approximate at least in the beginning the overlap from the data
%             if((true_corr_flag == TWO_SAMPLED) & (nsamples-3<=Nsamples/2)) % For now skip this !!!! 
%                 % Now do also sampling FROM THE DATA and compare to analytic result
%                 kept_frac_dist = sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, nsamples+d_var_n, Ngenes, alpha, iters, R, TRUE,miu,prior);
%                 samp_frac_mean_from_data(index, i) = mean(kept_frac_dist)/(alpha*Ngenes); samp_frac_std_from_data(index, i) = std(kept_frac_dist./(alpha*Ngenes));
%             end
%%%%            samp_frac_mean_from_data = samp_frac_mean; samp_frac_std_from_data = samp_frac_std; % Dummy wrong !!!! 

%             Now do also sampling and compare to analytic result
            i=i+1;

        end  % loop on nsamples_vec
        nsamples=max_nsamples;
        nsam_eff=1/(a*(nsamples+d_var_n)^b);
        sig = 1 / (Fisher_Zs_std*sqrt(nsam_eff));
%             max_y=100*(1+sqrt(sig));
%             dy=max_y/10000;
%             y=[0:dy:max_y];
        y=abs(Fisher_Zs/Fisher_Zs_std);
        P_in=(1 - normcdf( (saved_x_alpha-y)./sig )+normcdf( (-saved_x_alpha-y)./sig )).^2;

        
        % Now compute the whole distribution for the nsamples : 
        % Now show as an example the bahaviour for one nsamples value
        % This is called only twice : 
        if(true_corr_flag == TRUE_AND_SAMPLED) % i.e. true and sampled ..
            iters = 300; nsamples = 100; alpha  = alpha_vec(1); f_res = 0.025; soft_constrain_flag = 0;
            alpha
%             CheckPermTopFunc(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, nsamples-3, Ngenes, alpha, iters, f_res, R, soft_constrain_flag,miu,prior);
        end

        % adjust for Ngenes
        inf_limit_std(index,:) = inf_limit_std(index,:) ./ sqrt(Ngenes);


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



        figure; hold on; % Make new figure not to distroy the old one !!!! 
       
       % subplot(2,2,true_corr_flag+1); hold on; 

        % Analytic saddle-point approximation , here 'easy' stds.
  %%%      plot(nsamples_vec, inf_limit_frac(index,:) - stupid_independent_inf_limit_std(index,:), 'g+');
  %%%      plot(nsamples_vec, inf_limit_frac(index,:) + stupid_independent_inf_limit_std(index,:), 'g+');

        subplot(2,2,true_corr_flag+1); hold on;
        errorbar(nsamples_vec, samp_frac_mean(index,:), samp_frac_std(index,:), 'r');
        hold on,errorbar(nsamples_vec, data_samp_frac_mean(index,:), data_samp_frac_std(index,:), 'r');
        if(true_corr_flag == TRUE_AND_SAMPLED)
            errorbar(nsamples_vec, inf_limit_frac(index,:), inf_limit_std(index,:));
            legend('sampled from model','Saddle approx.', 2);
        else
            errorbar(nsamples_vec, samp_frac_mean_from_data(index,:), samp_frac_std_from_data(index,:), 'g');
            plot(nsamples_vec, inf_limit_frac(index,:));
            legend('sampled from model', 'sampled from data', 'Saddle approx.', 2);
        end


        plot_ind = plot_ind+1;

        index = index + 1;

        %%%      end  % loop on alpha
        %%% end  % loop on two samps flag


    end  % loop on true corrs


    subplot(2,2,3); hold on;  num_bins = max(250, floor(Ngenes/250)); % % num. of bins is arbitrary. Should be determined more cleverly
    [H bin_pos] = hist(SAVE_FISHER_Z{data_flag}, num_bins);
    hist(SAVE_FISHER_Z{data_flag}, num_bins); xlabel('Correlation'); ylabel('Freq.'); title([data_str ' Data Hist. of all genes Fisher Zs']);
    plot(Fisher_Zs_possible_vec, SAVE_FISHER_NORMAL_FIT_VEC{data_flag} .* Ngenes .* (bin_pos(2)-bin_pos(1)), 'r');
   
    % Set now wider axis
    AX = gca;    
    X_ax_lim = get(AX, 'xlim'); Y_ax_lim = get(AX, 'ylim');
    X_ax_lim(1) = min(X_ax_lim(1), X_ax_lim(1)*1.2);     X_ax_lim(2) = max(X_ax_lim(2), X_ax_lim(2)*1.2);
    Y_ax_lim(1) = min(Y_ax_lim(1), Y_ax_lim(1)*1.2);     Y_ax_lim(2) = max(Y_ax_lim(2), Y_ax_lim(2)*1.2);
    set(AX, 'xlim', X_ax_lim); set(AX, 'ylim', Y_ax_lim);
    
    
    data_mean = mean(SAVE_FISHER_Z{data_flag});
    data_std = std(SAVE_FISHER_Z{data_flag});
    saved_mean = SAVE_FISHER_Z_MEAN{data_flag};
    saved_std = SAVE_FISHER_Z_STD{data_flag};

    % Now take a random sample from the same Gaussian distribution and see
    % if it looks the same
    rand_gauss_samp = SAVE_FISHER_Z_STD{data_flag} .* randn(1, Ngenes) + SAVE_FISHER_Z_MEAN{data_flag};
    subplot(2,2,4); hold on;  num_bins = max(250, floor(Ngenes/250)); % % num. of bins is arbitrary. Should be determined more cleverly
    [H bin_pos] = hist(rand_gauss_samp, num_bins);
    hist(rand_gauss_samp, num_bins); xlabel('Correlation'); ylabel('Freq.'); title([data_str ' Data Hist. of rand. samp.']);
    plot(Fisher_Zs_possible_vec, SAVE_FISHER_NORMAL_FIT_VEC{data_flag} .* Ngenes .* (bin_pos(2)-bin_pos(1)), 'r');
    AX = gca;    
    set(AX, 'xlim', X_ax_lim); set(AX, 'ylim', Y_ax_lim);
    
    % Now do the plot that shows if indeed the Fisher Z's have std which is
    % independent of their value
    if( (calc_corrs_flag == TRUE) & (estimate_variance_flag == TRUE))
        corrs_iters = 2000;
        [N_VEC Pearson_corrs Fisher_Zs Rho_sig_mean Rho_sig_std Rho_bias_mean Rho_bias_std Z_sig_mean Z_sig_std Z_bias_mean ...
            Z_bias_std chunk_Rho_sig chunk_Z_sig chunk_Rho_bias chunk_Z_bias ] = ...
            CalcZVarianceFromDataFunc(R, alpha_vec, rand_flag, corrs_iters ); % take only max samples to save time
        subplot(2,3,4); hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. N=' num2str(N_VEC)]);
        xlabel('true P'); ylabel('std');
        subplot(2,3,5); hold on; plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. N='  num2str(N_VEC(1))]);
        xlabel('true Z'); ylabel('std');

        
        
        corrs_iters2=10;

        
        % Skip the following: Very time-consuming
        [N_VEC2 N_HALF_VEC2 Z_sig_mean_plot Q_sig_mean_plot ] = CalcTwoZVariancesFromDataFunc(R, alpha, rand_flag, corrs_iters2 );
               
         % Plot and compare both variances
         figure; hold on; plot(N_HALF_VEC2, Z_sig_mean_plot, '.'); plot(N_VEC2, Q_sig_mean_plot, '.r');  plot(N_HALF_VEC2, Q_sig_mean_plot(1:length(N_HALF_VEC2))-Z_sig_mean_plot, '.g'); 
         plot(N_HALF_VEC2, 1 ./ sqrt(N_HALF_VEC2-3), 'k');    plot(N_VEC2, std(Fisher_Zs) + 1 ./ sqrt(N_VEC2-3), 'k');
         legend('Ave Z std.', 'Q std.', 'diff', 'Fisher', 'Fisher plus one'); xlabel('Nsamples'); ylabel('Variance'); 
         title(['Std. of each genes Z and of total Q, Data is ' data_str]);

         
         
         % Now do fit
         
    end   % plot of std of the


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













