% Plot nicely all the gaussians q fits


TRUE = 1; FALSE = 0;

TOL = 0.00000001;


% data files
OLD_VANT_VEER = 1; NEW_ROSSETA = 2; WANG = 3;     % Breast gene expression
TOPIC = 4; % Author-topic matching
NIPS_DOROTHEA = 5; NIPS_ARCENE = 6; NIPS_DEXTER = 7; NIPS_GISETTE = 8; % Various datasetes from NIPS contest
ROSEN = 9; HEPATO =10; % New data's from Assif
Bhattacharjee=11;  BEER=12; % New Lung data
YEOH = 13; % leukemia
BRAIN = 14; KIM = 15; % New aging related datasets

RAND_DATA = 16; % Here we randomize a data so that we have control on it, and also can work with matlab 6.5 (no loading from files)

global BOOTSTRAP NON_OVERLAP MOG_NON_OVERLAP;
BOOTSTRAP=0; NON_OVERLAP=1; MOG_NON_OVERLAP=2;
sampling_flag = MOG_NON_OVERLAP;

global CONST_ALPHA LINEAR_ALPHA;
CONST_ALPHA = 0; LINEAR_ALPHA=1;
alpha_type = CONST_ALPHA; %%% LINEAR_ALPHA; % Choose whether alpha is kept constant or alpha*N_g ~ n

check_eps_delta_flag = 1;

% Choose the probability distribution of the TRUE corrleations
global UNIFORM GAUSSIAN LINEAR FROM_DATA mix_GAUSSIAN student_t;
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; FROM_DATA = 3; mix_GAUSSIAN=4;
student_t=5;% in the last one we take q simply according to the data bins
rand_flag = mix_GAUSSIAN;



% Choose if to take both top and bottom genes or just top
global ONE_SIDE TWO_SIDES;
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
global TRUE_AND_SAMPLED TWO_SAMPLED;
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TWO_SAMPLED;




kept_fraction=0.01;
KNN = 10; % parameter for knn missing values algorithm
Fisher_R_to_Z = 1; % Flag saying if to perform fisher RtoZ transformation
remove_zero_flag = TRUE;   % Flag saying to remove all the features with correlation
% zero from the data, since they cause bias.


load('all_datas_vars.mat');

data_flag = OLD_VANT_VEER; % NEW_ROSSETA; % NEW_ROSSETA;  % Choose which data to plot

calc_corrs_flag = TRUE; % calculate the correlation each time !!!!

index = 1; % Data index

figure; hold on; % Make one figure for all datasets

for data_flag =  [NEW_ROSSETA WANG ROSEN HEPATO Bhattacharjee BEER]
    % Here load the data
    if(data_flag == OLD_VANT_VEER)
        rand_flag = GAUSSIAN; % should be gaussian
        num_of_Gaussians = 2;
        load('../data/breastData.mat'); % The name of the data file containing the expression and labels.
        data_str = 'Breast Cancer [?]';
    end
    if(data_flag == NEW_ROSSETA)
        loading_big_rosseta = 1;
        load('../data/Rosetta_data.mat'); % The name of the data file containing the expression and labels.
        data_str = 'Breast Cancer [9]';

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
        data_str = 'Breast Cancer [10]';

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



    if(data_flag == ROSEN)
        load('../data/ROSEN.mat');
        data_str = 'Acute Lymphocitic Leukemia [4]';
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
        data_str = 'Lung Cancer [2]';
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
        a=c;
        b=d;
    end


    if(data_flag == HEPATO)
        load('../data/hepatocellular_carcinoma.mat');
        data_str = 'Hepatocellular Carcinoma [21]';
        rand_flag = GAUSSIAN;


        % % % %         % Liat's pre-processing
        % % % %         data(data < 30) = 30;
        % % % %         data = data + 0.001*rand(size(data,1),size(data,2)); %
        % % % %         data=log2(data);
        % % % %         present=real(pval<0.05);
        % % % %         present_idx=find(sum(present,2)>=size(pval,2)*kept_fraction);
        % % % %         present_data=data(present_idx,:);
        % % % %         R.dat=present_data;
        % % % %         R.Labels=labels';



        % My pre-processing
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



    if(data_flag == Bhattacharjee)
        rand_flag = GAUSSIAN;
        load('../data/Bhattacharjee_66.mat'); % The name of the data file containing the expression and labels.
        data_str = 'Lung Cancer [6]';
        %%%data(data<30)=30; % Lower the threshold
        data = data + 0*0.00001*randn(size(data,1),size(data,2));
        data=log2(data+1); % Try  ....
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

    curr_data_is = data_str

    Ngenes = size(R.dat, 1)
    Nsamples = size(R.dat, 2)



    % Compute the Pearsion correlation coefficients:
    if(calc_corrs_flag == TRUE)
        % Calculate the mean correlation of each gene with survival
        normed_labels = R.Labels - mean(R.Labels); normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
        normed_data = R.dat;
        %        normed_data(4965,:) = rand(1,Nsamples); % BADD THING !!!
        normed_data = normed_data - repmat(mean(normed_data, 2), 1, Nsamples);
        normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, Nsamples);
        %normed_data = normed_data';

        size(normed_data)
        size(normed_labels)

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


    % Do the hyperbolic tangent transform:
    Fisher_Zs = 0.5 * (log(1+Pearson_corrs) - log(1-Pearson_corrs));

    Fisher_Zs_mean = mean(Fisher_Zs); Fisher_Zs_std = std(Fisher_Zs);
    Fisher_Zs_min = min(Fisher_Zs);  Fisher_Zs_max = max(Fisher_Zs); Fisher_Zs_gap = Fisher_Zs_max-Fisher_Zs_min;
    Fisher_Zs_possible_vec = [(Fisher_Zs_min - 0.05*Fisher_Zs_gap):(Fisher_Zs_gap*0.01):(Fisher_Zs_max+0.05*Fisher_Zs_gap)];
    
    Fisher_Normal_fit_vec = (1/(sqrt(2*pi)*Fisher_Zs_std)) .* ...
        exp(-(Fisher_Zs_possible_vec-Fisher_Zs_mean).^2 ./ (2 .* Fisher_Zs_std^2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot the q_z distribution
    subplot(2,3,index);
    hold on;  num_bins = max(250, floor(Ngenes/250)); % % num. of bins is arbitrary. Should be determined more cleverly
    [H bin_pos] = hist(Fisher_Zs, num_bins);
    hist(Fisher_Zs, num_bins); xlabel('Z'); ylabel('P(Z)');  title(data_str);
    plot(Fisher_Zs_possible_vec, Fisher_Normal_fit_vec .* Ngenes .* (bin_pos(2)-bin_pos(1)), 'r');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y_min = 0; y_max = 220; x_min=-0.6; x_max = 0.6; AXIS([x_min x_max y_min y_max]); 

    
    index=index+1;
    
    clear R; clear dat; clear data;

end







% removed old crap