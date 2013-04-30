% This script assumes that we already learned an HMM model for each
% chromosome, and we now wish to compare them and build a distance matrix
% etc. etc. 


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

num_chroms = 22; 

do_log = 0; % Flag saying if to perform log transform
do_center_norm = 1; % Flag saying if to do centering and normalization

do_load = 0; % load from file to see what is best . 0 is risky !!! Only for the first time !!! 

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


max_iters = 100;  % maximum number of EM iterations
num_starts = 10;  % number of different starting points for EM
tolerance = 0.000001;  % tolerance for score improving 

seq_len = 1000; 
num_iters = 10; 

cv = 0;   % cross-validation technique - this is the number of samples we learn on each time. 
          % 0 means no cv - training and testing on same data (which is wrong)



for trials = 1:meta_run  % meta-run
    mean_x = []; % init to zero 
    std_x = []; % init to zero
    
    num_genes_in_chrom = [];  % Number of genes in each Chromosome 

    
    PI = {}; M = {}; N = {}; MEW = {}; SIGMA = {}; 
    
    for i=1:num_chroms    % First over all the chromosomes and load them .. 
        %%%%%%%%%%%%%   Load model i           
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
        load(cur_str);
        cd ..;
        
        % Copy parameters of the chromosoems 
        PI{i} = HMM_CHROM{2};
        M{i} = HMM_CHROM{3};
        N{i} = HMM_CHROM{4};
        MEW{i} = HMM_CHROM{5};
        SIGMA{i} = HMM_CHROM{6};            
        loaded_score_first_is = HMM_CHROM{7};
        %%%%%%%%%%%%%             
        
    end
    
    KL_DIST = zeros(num_chroms, num_chroms); 
    
    % Go over all couples to evaluate distances
    for i=1:num_chroms    % Go over all the chromosomes   
        now_doing_chrom = i
        for j=1:num_chroms  % Go over all the other chromosomes
            % Now do the comparison of the hmm models .... %%%%%%%%%%%%% Now call the HMM comparison module 
            KL_DIST(i, j) = ComputeHMMKLDistance(PI{i}, M{i}, N{i}, MEW{i}, SIGMA{i}, ...
                                    PI{j}, M{j}, N{j}, MEW{j}, SIGMA{j}, seq_len, num_iters);
        end
    end
                       
            
    figure; hold on; imagesc(KL_DIST); colorbar; title('Chromosomes Distance Matrix Between Models'); xlabel('chrom. num'); ylabel('chrom. num.');
    
    
    KL_DIST_SIM = (KL_DIST+KL_DIST')/2;
    figure; hold on; imagesc(KL_DIST_SIM); colorbar; title('Chromosomes SYMMETRIC Distance Matrix Between Models'); xlabel('chrom. num'); ylabel('chrom. num.');
    
end % meta-run




cd ../../