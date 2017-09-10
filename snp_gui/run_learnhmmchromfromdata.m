% Run the function for learning the models from the data. Here is what the
% Here set all the parameters
%
DataTypesVec  = {'Colon', 'Uveal', 'Downing', 'GBM', 'LEUKEMIA_SNPS', 'COLON_SNPS'};
COLON = 1;
UVEAL = 2;
DOWNING = 3;
GBM = 4;
LEUKEMIA_SNPS = 5;
COLON_SNPS = 6;
DataIndex = LEUKEMIA_SNPS;  % DOWNING COLON UVEAL

% Various data types
MRNA_EXP = 0; SNP_CHIPS=1;
DataType = MRNA_EXP;

% Downing Leukimia Data
if(DataIndex == DOWNING)
    full_path_data_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\downing\downing_U133_data_for_HMM.mat';
    full_path_model_and_output_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\downing\output\NewDowningHMMLogDifferentSizes.mat';
end

% Colon Cancer Data
if(DataIndex == COLON)
    full_path_data_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\colon\colon_data_for_HMM.mat';
    full_path_model_and_output_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\colon\output\NewColonHMMLogDifferentSizes.mat';
end

% Uveal Melanoma Data
if(DataIndex == UVEAL)
    full_path_data_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\uveal\Uveal_data_for_HMM.mat';
    full_path_model_and_output_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\uveal\output\NewUvealHMMLogDifferentSizes.mat';
end

% Glioblastoma Data
if(DataIndex == GBM)
    full_path_data_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\GBM\GBM_data_for_HMM.mat';
    full_path_model_and_output_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\GBM\output\NewGBMHMMLogDifferentSizes.mat';
end


% Here we load the new SNPs chips data:
if(DataIndex == LEUKEMIA_SNPS)

    sample_name = 'HD78_9_n';
    sample_name2 = 'HD78_9_d'; %%%% 'TEL74_5_n';

    DataType = SNP_CHIPS;
    chip_type = 'Hind'; genome_assembly = 'hg17';
    ChipStr = [chip_type '_annot_data_' genome_assembly '.mat'];
    full_path_data_file_name = ['..\data\Leukemia\' sample_name '_' chip_type];
    full_path_model_and_output_file_name = 'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\Leukemia\output\ALL_TEL_HD_SNP_HMM_LogDifferentSizes.mat';
    user_dir = '..\data\Leukemia';
    HMM_x_dim=3;
    SNPChipAnnotStruct = load(['..\Database\' ChipStr]); %, 'snp_ids', 'rs_ids', 'chr_loc_vec', 'chr_vec', 'strand');
    LDStruct = load(['..\Database\LD_' chip_type]); % load LD
    LD = LDStruct.LD;
    snp_ids = SNPChipAnnotStruct.snp_ids;
    rs_ids = SNPChipAnnotStruct.rs_ids;
    chr_loc_vec = SNPChipAnnotStruct.chr_loc_vec;
    chr_vec = SNPChipAnnotStruct.chr_vec;
    strand = SNPChipAnnotStruct.strand;
    SampleNames = {'HD78_9_n' 'HD78_9_d', 'HD82_3_n', 'HD82_3_d', 'HD86_7_n', 'HD86_7_d', 'HD90_1_n', 'HD90_1_d', ...
        'HD92_3_n', 'HD92_3_d', 'TEL74_5_n', 'TEL74_5_d', 'TEL76_7_n', 'TEL76_7_d', 'TEL80_1_n', 'TEL80_1_d', ...
        'TEL84_5_n', 'TEL84_5_d', 'TEL88_9_n', 'TEL88_9_d', 'TEL94_5_n', 'TEL94_5_d', 'TEL96_7_n', 'TEL96_7_d'  };
    SampleNames = {}; SampleNames{1} = 'TEL94_5_n';
    
    
%     SampleNames{1} = sample_name;
%     SampleNames{2} = sample_name2;
end
if(DataIndex == COLON_SNPS)

end




% home
%full_path_data_file_name = 'C:\Weizmann\hmm_chrom\hmm_chrom\data\colon\colon_data_for_HMM.mat';
%full_path_model_and_output_file_name = 'C:\Weizmann\hmm_chrom\hmm_chrom\data\colon\output\NewColonHMMLogDifferentSizes.mat';

num_chroms = 22;


HMM_x_dim = 3;  % Number of amplification levels
HMM_y_dim = 1;  % number of Gaussians for expression for each amplification level
learn_model_EM = 1;     % Flag saying if we want to learn the model
find_path_Viterbi = 1;  % Flag saying if we want to get the best pathes suggestions for amplifications
compute_chromosomes_distance_matrix = 0;   % Flag saying if to compute the distance matrix between chromosomes
use_locations = 0;     % Flag saying if to use distances between genes. Currently not recommended !
do_center_norm = 0;   % Flag saying if to do centering and normalization on data
do_fold_change = 0;   % Flag saying if to do fold-change with respect to some reference set of samples. Note ! Not recommended to do both fold-change and normalize !
do_log = 1;           % Flag saying if to do log. This should be done when using fold-change !
do_clip = 0;          % Clip the gaussian data into two levels (0, 1)
do_median_flag = 1;   % 1 - do median, 0 do mean (where do we do median/mean???)


do_determine_x_dim = 0;   % Flag saying if to determine X dimension automatically. In this case, the Maximal possible dimenstion
% is HMM_x_dim. The actual dimension that we will get is anything between 1 and HMM_x_dim



num_EM_iters = 50;          % Number of EM iterations. More iterations --> Better solution but increased running time
num_EM_starting_points = 10; % Number of EM starting points. More points --> Better solution but increased running time
EM_tolerance = 0.000001;    % tolerance to tell the EM to stop if likelihood improvements are too minor
meta_run = 1;               % How many times should we run the EM learning algorithm

num_KL_dist_iters = 10;      % Number of iterations when computing Kullback-Leibler distance between chromosomes
KL_dist_seq_len = 1000;     % Length of sequence to simulate when computing Kullback-Leibler distance between chromosomes

ChromosomesToRun = [1:22];   % Which chromosomes do you want to run on. Default is 1 to 22


run_flag = 1;   % flag saying if you want to run calculations
plot_flag = 0;  % flag saying if you want to plot results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% The plotting part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
ChromosomesToPlot = 21; %% ChromosomesToRun; %%% [21];   % Which chromosomes do you want to plot
PatientsToPlot = 1; %%% [1:22]; % [1 2 12];    % Which patients do you want plot

ExpressionPlot = 1;  % Plot the expression level
ViterbiPlot = 1;     % plot the most likely path
MarginalPlot = 1;    % plot the prob. for each amp. level

ChromosomesDistancePlot = 1; % Flag saying if to plot the chromosomal distnace matrix
num_samples_in_row = 4;      % Number of plots to show in every row
num_samples_in_column = 3;   % Number of plots to show in every column


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of User Controlls. Do not change anything below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%do_load = 1; % load from file to see what is best . 0 is risky since it ruins the already stored models !!! Use 0 Only for the first time !!!
%%% Note : In the first time we MUST run on all chromosomes !!!
if(exist(full_path_model_and_output_file_name))
    do_load = 1;
else
    do_load = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
HMMParamsStruct.x_dim = 3; HMMParamsStruct.y_dim = 1;
HMMParamsStruct.num_EM_iters = num_EM_iters;
HMMParamsStruct.num_EM_starting_points = num_EM_starting_points;
HMMParamsStruct.EM_tolerance = EM_tolerance;
HMMParamsStruct.learn_model_EM = learn_model_EM;
HMMParamsStruct.ChromosomesToRun = ChromosomesToRun;
HMMParamsStruct.ChromosomesToPlot = ChromosomesToRun; % Plot
HMMParamsStruct.derich = derich;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% The running part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if run_flag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Set a sample model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine parameters ...
    derich = 1.0/120.0; 
    P_chrom_copy_number_change = 0.001;  % Probability to change the copy number
    HMM_MODEL.SIGMA = [0.0778    0.2868    0.4103    0.3611    0.6931]; % From MOG fitting
    %%%%%%  0.645 * ones(1,2*HMMParamsStruct.x_dim-1); % Increase sigma ... remember that we have more than x_dim
    HMM_MODEL.N = ones(2*HMMParamsStruct.x_dim-1,1);
    HMM_MODEL.PI = ones(HMMParamsStruct.x_dim,1) ./ HMMParamsStruct.x_dim;
    HMM_MODEL.M = eye(HMMParamsStruct.x_dim) .* (1-2*P_chrom_copy_number_change);
    for j=1:HMMParamsStruct.x_dim-1
        HMM_MODEL.M(j+1,j) = P_chrom_copy_number_change;
        HMM_MODEL.M(j,j+1) = P_chrom_copy_number_change;
    end
    HMM_MODEL.M(1,1) = 1-P_chrom_copy_number_change;
    HMM_MODEL.M(HMMParamsStruct.x_dim,HMMParamsStruct.x_dim) = 1-P_chrom_copy_number_change;
    HMM_MODEL.MEW = [ 0.1422 0.7262 1.5627 2.3612 2.4599]'; % From MOG fitting %%%%[0.2 0.8 1.8 2.45 3]';  % Note: Mew here is of size 2*dim-1 !!!
    HMM_MODEL.PLACE_FLAG=1;
    HMM_MODEL.SPECIAL_MODELS_FLAG=1;
    
    
    
    HMM_MODEL = LDToPlaceM(LDStruct, HMM_MODEL, HMMParamsStruct); % update the PLACE_M matrix according to Linkage-Disequilibrium

    %%%%HMM_MODEL.PLACE_M = HMM_MODEL.PLACE_M*0 + 0.5;
    HMMParamsStruct.MODEL = HMM_MODEL; % Copy the model itself into the parameter's struct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End setting a sample model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for iter=1:meta_run
        % Now run the function
        [OutputFilesNames] = ...
            LearnHMMChromFromData(user_dir, SampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% The plotting part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if plot_flag
    for SinglePatientToPlot = PatientsToPlot
       DUMMY = PlotPatientOnChromResults(user_dir, SampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct);
%        Old Version - BADD        
%        DUMMY = PlotPatientOnChromResults(full_path_data_file_name, full_path_model_and_output_file_name, DataType, ...
%            ChromosomesToPlot, SinglePatientToPlot, ...
%            ExpressionPlot, ViterbiPlot, 1-MarginalPlot, do_fold_change, do_clip);
    end
end

