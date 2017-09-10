% Run the function for learning the models from the data. Here is what the
% Here set all the parameters
full_path_data_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\colon\colon_data_for_HMM.mat';
full_path_model_and_output_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\colon\output\colon_models_and_results_for_HMM_fold_change_on_modelsize_is_3.mat';

num_chroms = 22; 


HMM_x_dim = 3;  % Number of amplification levels
HMM_y_dim = 1;  % number of Gaussians for expression for each amplification level
learn_model_EM = 1;     % Flag saying if we want to learn the model
find_path_Viterbi = 1;  % Flag saying if we want to get the best pathes suggestions for amplifications
compute_chromosomes_distance_matrix = 1;   % Flag saying if to compute the distance matrix between chromosomes
use_locations = 0;     % Flag saying if to use distances between genes. Currently not recommended !
do_center_norm = 1;   % Flag saying if to do centering and normalization on data
do_fold_change = 0;   % Flag saying if to do fold-change with respect to some reference set of samples. Note ! Not recommended to do both fold-change and normalize !



num_EM_iters = 10;          % Number of EM iterations. More iterations --> Better solution but increased running time
num_EM_starting_points = 1; % Number of EM starting points. More points --> Better solution but increased running time
EM_tolerance = 0.000001;    % tolerance to tell the EM to stop if likelihood improvements are too minor
meta_run = 1;               % How many times should we run the EM learning algorithm

num_KL_dist_iters = 10;      % Number of iterations when computing Kullback-Leibler distance between chromosomes
KL_dist_seq_len = 1000;     % Length of sequence to simulate when computing Kullback-Leibler distance between chromosomes



do_load = 1; % load from file to see what is best . 0 is risky since it ruins the already stored models !!! Use 0 Only for the first time !!! 


run_flag = 1;   % flag saying if you want to run calculations
plot_flag = 1;  % flag saying if you want to plot results


% loading the needed data 
load(full_path_data_file_name);


% Here's what we got : 
% HMM_chromosome_arm   - which arm (e.g. 7q) every gene lies on
% HMM_data             - samples (expression data matrix ) of all genes
% HMM_genes            - labels of all genes
% HMM_locations        - locations (in nucleotides from chromosome starts of each gene on chromosome
% HMM_samples          - labels of all samples 
% HMM_ref_labels       - labels saying if this sample is in the reference ('Normal') group or in the set we want to check ('Cancer')




% % % b)  The first six letters are the Tissue Label, which identifies a unique microdissection of a specific frozen tissue block by adding a one letter Tissue Code.  We are excluding the letters I, L, O due to their potential to be mistaken for the numbers 1 and 0.  The key follows:
% % % Primary:
% % % A-- primary colon cancer
% % % B-- primary polyp
% % % C--primary polyp, high grade
% % % D--microadenoma
% % % E, F, G--other primary from same pt (see database for details)
% % % 
% % % Normal:
% % % H--normal mucosa
% % % J--normal lymph node
% % % K-- normal liver
% % % M--normal lung
% % % N, P, Q--other normal from same pt  (see database for details)
% % %  
% % % Metastasis:
% % % T--lymph node met
% % % U--liver met
% % % V--lung met
% % % W,X,Y--other met from same pt(see database for details)
% % % 


% Now compute the labeling ! 
NormalSamplesLabels = 'HJKMNPQ';
PrimaryTumorsSamplesLabels = 'ABCDEFG';
MetastasisSamplesLabels = 'TUVWXY';




HMM_ref_new_labels = zeros(1,length(HMM_ref_labels));

% do it in a primitive way ..
for i=1:length(HMM_ref_labels)
   
   if(isempty(find(NormalSamplesLabels == HMM_samples{i}(7))))
       HMM_ref_new_labels(i) = 1;
   end 
end

HMM_ref_labels = HMM_ref_new_labels;

save(full_path_data_file_name, 'HMM_chromosome_arm' , 'HMM_data', 'HMM_genes', 'HMM_locations' , 'HMM_samples', 'HMM_ref_labels');
