% % % Patients : 
% % % aAaM18408_LOH_6proa_8qai
% % % aBaM18383_LOH_6proa_8qai
% % % aCaM18037_LOH_6pnd_8qai
% % % aJM18664_LOH_6proa_8qroa
% % % aQaM19206_LOH_6proa_8qroa
% % % aRaM18670_LOH_6proa_8qai
% % % aSaM19219_LOH_6proa_8qai
% % % aTaM19227_LOH_6proa_8qai
% % % aUaM19210_LOH_6qai_8qai
% % % aVaM18402_LOH_6proa_8qai
% % % bDaM17570_NLOH_6pnd_8qai
% % % bFaM18017_NLOH_6pnd_8qai
% % % bGaM18385_NLOH_6proa_8qroa
% % % bIaM19241_NLOH_6qai_8qroa
% % % bJaM17401_NLOH_6pnd_8qai
% % % bLaM18672_NLOH_6proa_8qai
% % % bMaM19208_NLOH_6qai_8qroa
% % % bNaM19233_NLOH_6proa_8qroa
% % % bOaM18666_NLOH_6qai_8qroa
% % % bPaM19220_NLOH_6proa_8qroa

patient_names = {'aAaM18408_LOH_6proa_8qai',
    'aBaM18383_LOH_6proa_8qai',
    'aCaM18037_LOH_6pnd_8qai',
    'aJM18664_LOH_6proa_8qroa',
    'aQaM19206_LOH_6proa_8qroa',
    'aRaM18670_LOH_6proa_8qai',
    'aSaM19219_LOH_6proa_8qai',
    'aTaM19227_LOH_6proa_8qai',
    'aUaM19210_LOH_6qai_8qai',
    'aVaM18402_LOH_6proa_8qai',
    'bDaM17570_NLOH_6pnd_8qai',
    'bFaM18017_NLOH_6pnd_8qai',
    'bGaM18385_NLOH_6proa_8qroa',
    'bIaM19241_NLOH_6qai_8qroa',
    'bJaM17401_NLOH_6pnd_8qai',
    'bLaM18672_NLOH_6proa_8qai',
    'bMaM19208_NLOH_6qai_8qroa',
    'bNaM19233_NLOH_6proa_8qroa',
    'bOaM18666_NLOH_6qai_8qroa',
    'bPaM19220_NLOH_6proa_8qroa'};


% load data from some chromosome
path(path,'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom');  % weizmann
path(path,'C:\Weizmann\HMM_ido_2004_02_16\hmm_chromosome\hmm_chrom\hmm_chrom');  % home
path(path,'C:\Weizmann\HMM_ido_2004_02_16\hmm_chromosome\hmm_chrom');  % home

cd data\uveal;


plot_vec = 'byrgmck';  % Colors for plotting


%%% Temp because of no file : 


epsilon = 0.000000001;

num_chroms = 22; % 22;

do_log = 0; % Flag saying if to perform log transform
do_center_norm = 1; % Flag saying if to do centering and normalization

do_load = 1; % load from file to see what is best . 0 is risky !!! Only for the first time !!! 

place_flag = 0; % give a different distribution for every place. These are based on chromosome 3 
                % which is divided to 10 patients with LOH and 10 patients with NLOH

% load the start of the q-arm
% load('FirstLocationOnQArm.mat');


first_on_q =[142604331, 95203899, 94912966, 52697919, 50705877, 62387336, 63850117, 48697951,  66551451, 41965097, 55478369, ...
             39709262, 17055261, 18814688, 18962068, 46771070, 25766888, 17462263, 34390067, 30817204, 14665357, 14608115];
% for i=1:num_chroms
%     first_on_q(i) = 1; 
% end


meta_run = 1; % Number of runs to do 
x_dim = 2;  % two amplification levels
y_dim = 1;  % one gaussian for now

use_loc = 0; % flag saying if to use locations on chromosomes 


max_iters = 50;  % maximum number of EM iterations
num_starts = 5;  % number of different starting points for EM
tolerance = 0.000001;  % tolerance for score improving 

cv = 0;   % cross-validation technique - this is the number of samples we learn on each time. 
          % 0 means no cv - training and testing on same data (which is wrong)



for trials = 1:meta_run  % meta-run
    mean_x = []; % init to zero 
    std_x = []; % init to zero
    
    mean_chrom_vec = [];
    std_chrom_vec = [];
    
    
    num_genes_in_chrom = [];  % Number of genes in each Chromosome 
    for i=1:num_chroms    % Go over all the chromosomes 
        doing_chrom = i
        if(trials <= meta_run)
            close all;
        end
        
        % safety 
%         if(do_load == 0) 
%             do_load = 1;
%         end
%         
        HMM_CHROM = {}; % empty model 
        
        chr_cancer_exp = [];   % Start empty expression vec 
        chr_loc = [];
        
        chr_exp_arr = [];
        
        chr_str = ['_chr' num2str(i) '_'];
        chr_list = dir(['*' chr_str '*.mat']);    num_patients = length(chr_list);
        
        % Determine how many models are there 
        if(cv == 0)
            num_cv_models = 1;
            train_set_size = num_patients; test_set_size = num_patients;
        else
            num_cv_models = num_patients/cv;  % we assume that this is integer !!!
            train_set_size = num_patients - cv; test_set_size = cv;
        end
        
        
        chr_list_disp = chr_list; % used for display .. 
        
        % Now make good names ... 
        for(j = 1:length(chr_list)) % Go over to find a matching patient 
            cur_number = patient_names{j}(5:9);
            for(k = 1:length(chr_list)) % Go over to find a matching patient again
                if(~isempty(findstr(chr_list(k).name, cur_number)))
                    chr_list_disp(k).name = patient_names{j};
                end        
            end
        end
        
        for(j = 1:length(chr_list)) % Go over to find a matching patient  to avoid underscores
            chr_list_disp(j).name(findstr(chr_list_disp(j).name, '_')) = '-';
        end
        
        % Now permute the chr_list in order for the cross-validation to
        % take one from each class !!! 
        permy = reshape([1:10; 11:20], 1, 20); 
        permed_chr_list = chr_list(permy);
        permed_chr_list_disp = chr_list_disp(permy); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Here read once everything !!! 
        all_chr_exp_arr = [];
        all_chr_cancer_exp = [];
        all_chr_loc = [];
        for(j = 1:num_patients) % loop over patients for training !!! 
            load(permed_chr_list(j).name);  
            
            if(do_log)
                all_chr_exp_arr(j,:) = log(max(cancer_exp_vec', 1));  % take one by one with log !!! 
            else
                all_chr_exp_arr(j,:) = cancer_exp_vec';  % take one by one !!! 
            end
            chr_loc_cell{j} = loc_vec';   
            
            % Now add to everybody 
            all_chr_cancer_exp = [all_chr_cancer_exp cancer_exp_vec'];
            all_chr_loc = [all_chr_loc loc_vec'];
        end
        
        num_genes_in_chrom(i) = length(cancer_exp_vec);  % Set the number of genes
        
        if(do_log)
            all_chr_cancer_exp = log(max(all_chr_cancer_exp', 1));  % Do log transform        
        else
            all_chr_cancer_exp = all_chr_cancer_exp';  % No log transform
        end
        
        % take before normalizing 
        mean_chrom_vec(i) = mean(all_chr_cancer_exp);
        std_chrom_vec(i) = std(all_chr_cancer_exp);
        
        
        % Do centering and normalization 
        if(do_center_norm)    
            all_mean_vec = mean(all_chr_exp_arr);
            all_std_vec = std(all_chr_exp_arr);            
            all_std_vec = max(all_std_vec, epsilon); % avoid zero std. 
            
            all_chr_exp_arr = all_chr_exp_arr - repmat(all_mean_vec, size(all_chr_exp_arr, 1), 1);  % subtruct mean 
            all_chr_exp_arr = all_chr_exp_arr ./ repmat(all_std_vec, size(all_chr_exp_arr, 1), 1);  % div by standard deviation
            
            for j = 1:length(permed_chr_list)  % loop over patients only in the training set
                all_chr_cancer_exp(1+(j-1)*size(all_chr_exp_arr, 2):j*size(all_chr_exp_arr, 2)) = ...
                    (all_chr_cancer_exp(1+(j-1)*size(all_chr_exp_arr, 2):j*size(all_chr_exp_arr, 2)) - all_mean_vec') ./ all_std_vec';
            end
        end

        
        % Here done reading          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        for(cv_iter = 1:num_cv_models)   % run different models for cross-validation
                    
            % empty everything everytime
            chr_cancer_exp = [];   % Start empty expression vec 
            chr_loc = [];            
            chr_exp_arr = [];
            
            
            % Set the training and test sets 
            cv_train_set = setdiff([1:num_patients], [(cv_iter-1)*cv+1:cv_iter*cv]);    
            if(cv == 0)
                cv_test_set = cv_train_set;
            else
                cv_test_set = [(cv_iter-1)*cv+1:cv_iter*cv];  
            end
            
            chr_exp_arr = all_chr_exp_arr(cv_train_set,:);  % NEWWWWW !!!!! 
            
            for j = cv_train_set
                chr_loc = [chr_loc  all_chr_loc((j-1)*num_genes_in_chrom(i) + 1:j*num_genes_in_chrom(i))];
                chr_cancer_exp = [chr_cancer_exp  all_chr_cancer_exp((j-1)*num_genes_in_chrom(i) + 1:j*num_genes_in_chrom(i))'];                
            end
            
            % Do transpose 
            chr_loc = chr_loc';
            chr_cancer_exp = chr_cancer_exp'; 
            
            
            % save the originals !!! 
            orig_exp_arr = chr_exp_arr; 
            orig_cancer_exp = chr_cancer_exp; 
            
            % Do centering and normalization 
            if(do_center_norm)    
                mean_vec = mean(chr_exp_arr);
                std_vec = std(chr_exp_arr);            
                std_vec = max(std_vec, epsilon); % avoid zero std. 
                
                chr_exp_arr = chr_exp_arr - repmat(mean_vec, size(chr_exp_arr, 1), 1);  % subtruct mean 
                chr_exp_arr = chr_exp_arr ./ repmat(std_vec, size(chr_exp_arr, 1), 1);  % div by standard deviation
                
                jj=0;
                for j = cv_train_set %%1:length(permed_chr_list))  % loop over patients only in the training set
                    jj=jj+1; % inner index
                    size(chr_cancer_exp(1+(jj-1)*size(chr_exp_arr, 2):jj*size(chr_exp_arr, 2)));
                    size(mean_vec');
                    size(std_vec');
                    chr_cancer_exp(1+(jj-1)*size(chr_exp_arr, 2):jj*size(chr_exp_arr, 2)) = ...
                        (chr_cancer_exp(1+(jj-1)*size(chr_exp_arr, 2):jj*size(chr_exp_arr, 2)) - mean_vec') ./ std_vec';
                end
            end
            
            
            
            % Use place gaussians : compute their mean and variance
            if(place_flag)
                mean_vec = [mean(chr_exp_arr(1:train_set_size/2,:))' mean(chr_exp_arr(train_set_size/2+1:train_set_size,:))'];
                std_vec = [std(chr_exp_arr(1:train_set_size/2,:))' std(chr_exp_arr(train_set_size/2+1:train_set_size,:))'];            
                std_vec = max(std_vec, epsilon); % avoid zero std.                                        
            end
            
            
%             figure; hold on;
%             plot(mean(chr_exp_arr(1:train_set_size/2,:)), '.'); plot(mean(chr_exp_arr(train_set_size/2+1:train_set_size,:)), 'r.'); title(['mean of every gene chrom' i]); legend('LOH pat','NLOH pat');
            
            
            count_smaller = size(find(mean(chr_exp_arr(1:train_set_size/2,:)) > mean(chr_exp_arr(train_set_size/2+1:train_set_size,:))));
            count_bigger = size(find(mean(chr_exp_arr(1:train_set_size/2,:)) <= mean(chr_exp_arr(train_set_size/2+1:train_set_size,:))));
            
       %%     chr_loc = chr_loc';
            
            
            
            % Train the model 
            mean_vec_size = size(mean_vec);
            std_vec_size = size(std_vec);
            x_dim_isss = x_dim;
            y_dim_isss = y_dim;
            
            
            
            mean_vec_rep = repmat(mean_vec, num_patients, 1);
            std_vec_rep = repmat(std_vec, num_patients, 1);
            mean_size_rep = size(mean_vec_rep);
            std_size_rep = size(std_vec_rep);
            
            %   break;
            size(chr_cancer_exp)
            size(chr_loc)
            
            
            [PI M N MEW SIGMA LogScore] = ...
                TrainHMMFromDataEM(chr_cancer_exp, chr_loc, use_loc, x_dim, y_dim, ...
                place_flag, mean_vec_rep, std_vec_rep, ...
                max_iters, num_starts, tolerance);
            
            current_score_is = LogScore
            
            
            
            % Now compare to the best thing we have remembered !! 
            cur_str = [ 'hmm_model_chr_' num2str(i) ];
            if(place_flag)
                cur_str = [cur_str '_place'];
            end
            if(use_loc)
                cur_str = [cur_str '_use_loc'];
            end
            
            if(cv > 0)  % here cross-validation
                cur_str = [cur_str '_cv_out' num2str(cv_test_set(1)) '-' num2str(cv_test_set(end))];                  
            end
            
            cd new_output;   % Transfer to output directory 
            if(do_load)
                load(cur_str);
                
                if(LogScore > HMM_CHROM{7})  % Score has improved !!! 
                    
                    % Save the model 
                    HMM_CHROM{1} = i;  
                    HMM_CHROM{2} = PI;
                    HMM_CHROM{3} = M;
                    HMM_CHROM{4} = N;
                    HMM_CHROM{5} = MEW;
                    HMM_CHROM{6} = SIGMA;
                    HMM_CHROM{7} = LogScore;
                    
                    save(cur_str, 'HMM_CHROM')
                else            
                    PI = HMM_CHROM{2};
                    M = HMM_CHROM{3};
                    N = HMM_CHROM{4};
                    MEW = HMM_CHROM{5};
                    SIGMA = HMM_CHROM{6};            
                    loaded_score_is = HMM_CHROM{7}
                end    
                
            else
                % Save the model anyway. Be carefull !!! 
                HMM_CHROM{1} = i;  
                HMM_CHROM{2} = PI;
                HMM_CHROM{3} = M;
                HMM_CHROM{4} = N;
                HMM_CHROM{5} = MEW;
                HMM_CHROM{6} = SIGMA;
                HMM_CHROM{7} = LogScore;
                
                % Saving !!! 
                saving = cur_str 
                
                save(cur_str , 'HMM_CHROM')
            end 
            
            
            % Find Viterbi best path
            CalllingViter = 1111
            [Viterbi_Path Gamma_Probs] = FindBestPathViterbi(chr_cancer_exp, chr_loc, use_loc, PI, M, N, MEW, SIGMA, ...
                place_flag, mean_vec_rep, std_vec_rep);  % add the place variables 
            %%%%    [Viterbi_Path Gamma_Probs] = FindBestPathViterbi(log_exp_vec, loc_vec, PI, M, N, MEW, SIGMA); 
            
            % Find the true mean and standard deviation of the chromosome (before
            % normalization)
            
            for x = 1:x_dim % go over x levels
                mean_x(x,i) = sum(Gamma_Probs(x,:)' .* orig_cancer_exp) / sum(Gamma_Probs(x,:)); %%% length(cancer_exp_vec); 
                std_x(x,i) = sum(Gamma_Probs(x,:)' .* (orig_cancer_exp-mean_x(x,i)) .* (orig_cancer_exp-mean_x(x,i))) / sum(Gamma_Probs(x,:));
            end     
            
            should_be_zero = sum(sum(Gamma_Probs)) - length(chr_cancer_exp)
            
            cd ..;  % go back to input files
        end % of cv loop over different sets of patients 
        
        cd new_output;   % Go to output dir
        
        figure;   % Do one plot for all the patients
        subplot(4,5,1); title(['Chromosome : ' num2str(i)], 'fontsize', 5);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Now do some more loops on cv sets for output ...
        for(cv_iter = 1:num_cv_models)   % run different models for cross-validation
            
            % Set the training and test sets 
            cv_train_set = setdiff([1:num_patients], [(cv_iter-1)*cv+1:cv_iter*cv]);    
            if(cv == 0)
                cv_test_set = cv_train_set;
            else
                cv_test_set = [(cv_iter-1)*cv+1:cv_iter*cv];  
            end

            
            jj=0;
            for j=cv_test_set %%% j = 1:length(permed_chr_list)  % loop on patients again
                jj=jj+1;  % count them one by one
                subplot(4, 5, permy(j)); hold on; 
                
                AXIS([0 max(chr_loc_cell{j})*1.1 1.2*min(all_chr_exp_arr(j,:)) 1.2*max(all_chr_exp_arr(j,:))])    %    AXIS([XMIN XMAX YMIN YMAX])
                
                plot( chr_loc_cell{j}, all_chr_exp_arr(j,:), '.');  ylabel('Exp. Level', 'fontsize', 5);   
                
                H = line([first_on_q(i), first_on_q(i)], [1.2*min(all_chr_exp_arr(j,:)), 1.2*max(all_chr_exp_arr(j,:))]); % set the boundary between p and q arms
                set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
                if(j == 18)
                    xlabel(['Loc. Chrom. '  num2str(i)] , 'fontsize', 5);
                end
                title(['Pat. ' permed_chr_list_disp(j).name ], 'fontsize', 5);  %%% 'Chromosome ' num2str(i) ' with  ' num2str(x_dim) ' levels ' num2str(y_dim) ' mixtures'] );         
                
                
            end % loop on patients again
        end  % loop on cv again
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        print ('-djpeg'  ,['UvealExpressionChrom-' num2str(i)]); % print ('-depsc2'  ,['UvealExpressionChrom-' num2str(i)]);
        
        figure;   % Do one plot for all the patients
        subplot(4,5,1); title(['Chromosome : ' num2str(i)], 'fontsize', 5);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for(cv_iter = 1:num_cv_models)   % run different models for cross-validation
            
            % Set the training and test sets 
            cv_train_set = setdiff([1:num_patients], [(cv_iter-1)*cv+1:cv_iter*cv]);    
             if(cv == 0)
                cv_test_set = cv_train_set;
            else
                cv_test_set = [(cv_iter-1)*cv+1:cv_iter*cv];  
            end

            for j = cv_test_set  %%1:length(permed_chr_list))  % loop on patients again
% % %                 doing_now_patient = j
% % %                 size(all_chr_exp_arr(j,:))
% % %                 size(chr_loc_cell{j})
% % %                 size(mean_vec)
% % %                 size(std_vec)
                [Viterbi_Path_Patient Gamma_Probs_Patient] = FindBestPathViterbi(all_chr_exp_arr(j,:), chr_loc_cell{j}, use_loc, PI, M, N, MEW, SIGMA, ...
                    place_flag, mean_vec, std_vec);  % add the place variables all_chr_exp_arr(j,:); 
                %%%
                %%%     figure;
                subplot(4, 5, permy(j)); hold on; 
                AXIS([0 max(chr_loc_cell{j})*1.1 0 1])    %    AXIS([XMIN XMAX YMIN YMAX])
                
                %%%      subplot(2,1,1);
                plot( chr_loc_cell{j}, Viterbi_Path_Patient, '.');  ylabel('Amp. Level', 'fontsize', 5);   
                
                H = line([first_on_q(i), first_on_q(i)], [0, x_dim-1]); % set the boundary between p and q arms
                set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
                if(j == 18)
                    xlabel(['Loc. Chrom. '  num2str(i)], 'fontsize', 5 );
                end
                title(['Pat. ' permed_chr_list_disp(j).name ], 'fontsize', 5);  %%% 'Chromosome ' num2str(i) ' with  ' num2str(x_dim) ' levels ' num2str(y_dim) ' mixtures'] );                      
                
            end % loop on patients again
        end  % loop on cv again
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        print ('-djpeg'  ,['UvealViterbiChrom-' num2str(i)]);  %       print ('-depsc2'  ,['UvealViterbiChrom-' num2str(i)]);
        
        
        a=figure;   % Do one plot for all the patients, now with gammas
        subplot(4,5,1); title(['Chromosome : ' num2str(i)], 'fontsize', 5);
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for(cv_iter = 1:num_cv_models)   % run different models for cross-validation
            
            % Set the training and test sets 
            cv_train_set = setdiff([1:num_patients], [(cv_iter-1)*cv+1:cv_iter*cv]);    
             if(cv == 0)
                cv_test_set = cv_train_set;
            else
                cv_test_set = [(cv_iter-1)*cv+1:cv_iter*cv];  
            end

            for(j = 1:length(permed_chr_list))  % loop on patients again
                
                [Viterbi_Path_Patient Gamma_Probs_Patient] = FindBestPathViterbi(all_chr_exp_arr(j,:), chr_loc_cell{j}, use_loc, PI, M, N, MEW, SIGMA, ...
                    place_flag, mean_vec, std_vec);  % add the place variables ); 
                subplot(4, 5, permy(j)); hold on; 
                
                
                Gamma_Probs_Patient_Cum = cumsum(Gamma_Probs_Patient);                        
                
                
                
                AXIS([0 max(chr_loc_cell{j})*1.1 0 1.1])    %    AXIS([XMIN XMAX YMIN YMAX])
                hold on; 
                for k=1:x_dim
                    %    bar(Gamma_Probs_Patient_Cum(i,:), ['.' plot_vec(i) ]);    
                    plot(chr_loc_cell{j}, Gamma_Probs_Patient_Cum(k,:), ['.' plot_vec(k) ]);    
                end
                if(j == 18)  % display in the middle ! 
                    xlabel(['Loc. Chrom. '  num2str(i)], 'fontsize', 5 );
                end            
                % xlabel('Chromosome location'); 
                ylabel('Amp. Prob.', 'fontsize', 5); 
                title(['Pat. ' permed_chr_list_disp(j).name ], 'fontsize', 5); %%% [''Pat. 'Marginal Predicted Amplifications ' num2str(x_dim) ' levels ' num2str(y_dim) ' mixtures'] ); 
                H = line([first_on_q(i), first_on_q(i)], [0, x_dim-1+0.1]); % set the boundary between p and q arms
                set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
                
                
            end % loop on patients again
        end  % loop on cv again
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        print ('-djpeg'  ,['UvealMarginalChrom-' num2str(i)]);   %      print ('-depsc2'  ,['UvealMarginalChrom-' num2str(i)]);
        
        cd ..; % Go back to input files 
    end  % loop on chromosomes 
    
end % meta-run



% Now do some chromosomal plots 

mew_low_vec = []; mew_high_vec = [];

sigma_low_vec = []; sigma_high_vec = [];

high_to_low_vhr_vec = []; low_to_high_chr_vec = [];



cd new_output;   % Transfer to output directory again
for i=1:num_chroms
    cur_str = [ 'hmm_model_chr_' num2str(i) ];
    if(place_flag)
        cur_str = [cur_str '_place'];
    end
    if(use_loc)
        cur_str = [cur_str '_use_loc'];
    end
    
    if(cv > 0)  % here cross-validation
        cur_str = [cur_str '_cv_out' num2str(cv_test_set(1)) '-' num2str(cv_test_set(end))];                  
    end
    
    load(cur_str);
    mew_low_chr_vec(i) = HMM_CHROM{5}(1);
    mew_high_chr_vec(i) = HMM_CHROM{5}(2);
    
    sigma_low_chr_vec(i) = HMM_CHROM{6}(1);
    sigma_high_chr_vec(i) = HMM_CHROM{6}(2);
    
    high_to_low_vhr_vec(i) = HMM_CHROM{3}(1,2);
    low_to_high_vhr_vec(i) = HMM_CHROM{3}(2,1);
    
end
cd ..;   % transfer back again



figure; 
subplot(2,2,1);
hold on; plot(mew_low_chr_vec, '+'); plot(mew_high_chr_vec, '+r'); legend('low level', 'high level');
xlabel('Chromosome'); ylabel('Mean Expression Level (normalized)'); title('Mean Expression level for two amplification levels for each Chromosome - Normalized');

subplot(2,2,2);
hold on; plot(mean_x(1,:), '+'); plot(mean_x(2,:), '+r'); legend('low level', 'high level');
xlabel('Chromosome'); ylabel('Mean Expression Level (original)'); title('Mean Expression level for two amplification levels for each Chromosome - Original');

subplot(2,2,3);
hold on; plot(sigma_low_chr_vec, '+'); plot(sigma_high_chr_vec, '+r'); legend('low level', 'high level');
xlabel('Chromosome'); ylabel('Std. Deviation Expression Level (normalized)'); title('Std. Deviation Expression level for two amplification levels for each Chromosome - Normalized');

subplot(2,2,4);
hold on; plot(std_x(1,:), '+'); plot(std_x(2,:), '+r'); legend('low level', 'high level');
xlabel('Chromosome'); ylabel('Std. Deviation Expression Level (original)'); title('Std. Deviation Expression level for two amplification levels for each Chromosome - Original');



figure; hold on; plot(high_to_low_vhr_vec, '+'); plot(low_to_high_vhr_vec, '+r'); legend('high to low', 'low to high');
xlabel('Chromosome'); ylabel('Transfer Probabilities'); title('Prob. of transfer for each amp. level for each Chromosome');


% Now check the differences in the locations

loc_diff = chr_loc_cell{1}(2:end) - chr_loc_cell{1}(1:end-1);
figure; plot(loc_diff, '.');
figure; hist(loc_diff, 100);











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now compute distnaces between chromosomes
KL_chrom_dist = zeros(num_chroms, num_chroms); 
for i=1:num_chroms
    for j=1:num_chroms
        KL_chrom_dist(i, j) = log(std_chrom_vec(j)/std_chrom_vec(i))  + ...
            ((mean_chrom_vec(i)-mean_chrom_vec(j))^2 + std_chrom_vec(i)^2-std_chrom_vec(j)^2) / (2*std_chrom_vec(j)^2);
    end
end


KL_chrom_dist_sim = (KL_chrom_dist+KL_chrom_dist')/2;
figure; hold on; imagesc(KL_chrom_dist); colorbar; title('Chromosomes Distance Matrix Between Expression Distributions'); xlabel('chrom. num'); ylabel('chrom. num.');
figure; hold on; imagesc(KL_chrom_dist_sim); colorbar; title('Chromosomes SYMMETRIC Distance Matrix Between Expression Distributions'); xlabel('chrom. num'); ylabel('chrom. num.');



cd ../../