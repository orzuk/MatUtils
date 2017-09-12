% An interface with the C function SimulateSequenceFromModel. 
function [OUT_VEC HIDDEN_VEC] = SimulateSequenceFromModelMatlab( full_path_model_and_output_file_name, ChromosomeToPlot, KMin, num_points)


% loading the needed data 
%load(full_path_data_file_name);


% Here's what we got : 
% HMM_chromosome_arm   - which arm (e.g. 7q) every gene lies on
% HMM_data             - samples (expression data matrix ) of all genes
% HMM_genes            - labels of all genes
% HMM_locations        - locations (in nucleotides from chromosome starts of each gene on chromosome
% HMM_samples          - labels of all samples 
% HMM_ref_labels       - labels saying if this sample is in the reference ('Normal') group or in the set we want to check ('Cancer')


load(full_path_model_and_output_file_name);

% Here's what we got : 

% HMM_MODEL_SAVED     - The markov models
% Viterbi_Path        - The Viterbi best pathes for each Patient and each Chromosome
% Gamma_Probs         - The marginal probabilities for each Patient and each Chromosome 
% HMM_CHROM_KL_DIST   - The Chromosomal Distance Matrix

OUT_VEC = SimulateSequenceFromModel(HMM_MODEL_SAVED{ChromosomeToPlot}{KMin}.PI, HMM_MODEL_SAVED{ChromosomeToPlot}{KMin}.M, HMM_MODEL_SAVED{ChromosomeToPlot}{KMin}.N, ...
                          HMM_MODEL_SAVED{ChromosomeToPlot}{KMin}.MU, HMM_MODEL_SAVED{ChromosomeToPlot}{KMin}.SIGMA, num_points);



