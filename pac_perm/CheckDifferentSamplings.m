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
Bhattacharjee=11; YEOH = 12; % New Lung data
BRAIN = 13; KIM = 14; % New aging related datasets

RAND_DATA = 15; % Here we randomize a data so that we have control on it, and also can work with matlab 6.5 (no loading from files)


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
iters = 50;



% Now chose which data to load - New : Run on all data's !!!!

num_samples_needed = zeros(1,2);
std_needed = zeros(1,2);


SMALL_GENES_NUM = 1000; % Take only this amount of genes ...

Fisher_R_to_Z = 1; % Flag saying if to perform fisher R-to-Z transformation
one_side_flag = TWO_SIDES; % make sure we do two-sided unless later we change the flag

clear R; % free memory



data_iters = 50; % How many times to generate data!


f_overlap = zeros(1, data_iters); f_overlap_single = zeros(1, data_iters);



% First step is short : 
for data_iter = 1:data_iters

    data_iter_is = data_iter    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here generate R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DISCRETE_LABELS = 1; CONTINUOUS_LABELS = 0;% 1 discrete 0/1 0 Gaussians

    labels_type = CONTINUOUS_LABELS;

    Ngenes = 5000; Nsamples = 30;

    if(labels_type == DISCRETE_LABELS)
        R.Labels = round(rand(1, Nsamples)); % 0/1 labels
        pos_ind = find(R.Labels == 1);
    else
        R.Labels = randn(1, Nsamples); % Gaussian labels
    end

    rand_flag = GAUSSIAN;


    % Generate a more sophisticated  random model:
    Q_corr_vec = 1*(randn(1, Ngenes)); % Correlation of each gene with the survival
    W_corr_vec = 1*(randn(1, Ngenes));

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

    Ngenes = size(R.dat, 1)
    Nsamples = size(R.dat, 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here finished generating R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Check the width of Q in this specific experiment
    
    true_q_width_1st = std(Q_corr_vec)
    curr_q_width_1st = std(atanh(corr(R.dat', R.Labels')))
    
    q_width(data_iter) = std(atanh(corr(R.dat', R.Labels'))); 
    
    
    
    coeffs_vec = []; corr_mat = [];
    
    % Now we need to sample from the data : 
    [kept_frac_dist f_all_genes_from_data f_all_genes_from_data_one_NTOP] = ...
        sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, 1111, 15, Ngenes, alpha, iters, R, ...
           TRUE, miu, prior, 0.000000001 , coeffs_vec, corr_mat);
    
    f_overlap(data_iter) = mean(kept_frac_dist)/(alpha*Ngenes);  
           
    
    
    [kept_frac_dist f_all_genes_from_data f_all_genes_from_data_one_NTOP] = ...
        sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, 1111, 15, Ngenes, alpha, 1, R, ...
           TRUE, miu, prior, 0.000000001 , coeffs_vec, corr_mat);
           
    f_overlap_single(data_iter) = mean(kept_frac_dist)/(alpha*Ngenes);  


end % Loop on generating data's




q_ave_width = mean(q_width);
q_ave_width_ratio = q_ave_width / q_width(end); 


save_q_width = q_width; % Save for the next phase ...

% Now plot the results : 

figure; 
subplot(2,2,1); hold on; hist(f_overlap);  title(['non-over sets. mean: ' num2str(mean(f_overlap)) ' std: ' num2str(std(f_overlap))]); 
subplot(2,2,2); hist(f_overlap_single);  title(['f dist. mean: ' num2str(mean(f_overlap_single)) ' std: ' num2str(std(f_overlap_single))]); 
subplot(2,2,3); plot(f_overlap, f_overlap_single, '.'); title('one samp. vs. many non-overlapping samps.'); xlabel('many samps.'); ylabel('one samp');
subplot(2,2,4); plot(q_width, f_overlap, '.'); title('Dependence of overlap on width of Q'); xlabel('Q width'); ylabel('Overlap');



q_min = quantile(save_q_width, 0.4); q_max = quantile(save_q_width, 0.6);


% Now do the long part. Here we put our data in the center 
% First step is short : 

f_overlap = zeros(1, data_iters); f_overlap_single = zeros(1, data_iters);


for data_iter = 1:data_iters

    data_iter_2nd_time_is = data_iter   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here generate R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DISCRETE_LABELS = 1; CONTINUOUS_LABELS = 0;% 1 discrete 0/1 0 Gaussians

    labels_type = CONTINUOUS_LABELS;

    Ngenes = 5000; Nsamples = 30;



    rand_flag = GAUSSIAN;


    
    randomize_flag=0;
    
    while(randomize_flag == 0)
        
        if(labels_type == DISCRETE_LABELS)
            R.Labels = round(rand(1, Nsamples)); % 0/1 labels
            pos_ind = find(R.Labels == 1);
        else
            R.Labels = randn(1, Nsamples); % Gaussian labels
        end

        % Generate a more sophisticated  random model:
        Q_corr_vec = 1*(randn(1, Ngenes)); % Correlation of each gene with the survival
        W_corr_vec = 1*(randn(1, Ngenes));

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

        Ngenes = size(R.dat, 1)
        Nsamples = size(R.dat, 2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here finished generating R
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        true_q_width = std(Q_corr_vec)
        curr_q_width = std(atanh(corr(R.dat', R.Labels')))
        q_min
        q_max

        randomize_flag = 1;
        if(curr_q_width < q_min) 
           randomize_flag = 0;
        end
        if(curr_q_width > q_max) 
           randomize_flag = 0;
        end   
        
        % Check the width of Q in this specific experiment
        q_width(data_iter) = std(atanh(corr(R.dat', R.Labels')));


    end % while

    coeffs_vec = []; corr_mat = [];

    % Now we need to sample from the data :
    [kept_frac_dist f_all_genes_from_data f_all_genes_from_data_one_NTOP] = ...
        sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, 1111, 15, Ngenes, alpha, iters, R, ...
           TRUE, miu, prior, 0.000000001 , coeffs_vec, corr_mat);
    
    f_overlap(data_iter) = mean(kept_frac_dist)/(alpha*Ngenes);  
           
    
    
    [kept_frac_dist f_all_genes_from_data f_all_genes_from_data_one_NTOP] = ...
        sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, 1111, 15, Ngenes, alpha, 1, R, ...
           TRUE, miu, prior, 0.000000001 , coeffs_vec, corr_mat);
           
    f_overlap_single(data_iter) = mean(kept_frac_dist)/(alpha*Ngenes);  


end % Loop on generating data's



time_elapsed = cputime -ttt_time










