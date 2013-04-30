% Plot a mixture model and see how does it look !!!! 
function [COMP_HMM_MODELS] = CheckedClippedData( full_path_data_file_name, full_path_model_and_output_file_name)
load(full_path_data_file_name);


% Here's what we got : 
% HMM_chromosome_arm   - which arm (e.g. 7q) every gene lies on
% HMM_data             - samples (expression data matrix ) of all genes
% HMM_genes            - labels of all genes
% HMM_locations        - locations (in nucleotides from chromosome starts of each gene on chromosome
% HMM_samples          - labels of all samples 
% HMM_ref_labels       - labels saying if this sample is in the reference ('Normal') group or in the set we want to check ('Cancer')

load(full_path_model_and_output_file_name);


for i=ChromosomesToPlot
    ClippedChromData = ClipModelObservedData( HMM_MODEL_SAVED{i}{2}, HMM_chrom_data);
            
    PlotPatientOnChromResults(full_path_data_file_name, full_path_model_and_output_file_name, ChromosomesToPlot, PatientsToPlot, ...
                                      ExpressionPlot, ViterbiPlot, 1-MarginalPlot, do_fold_change);
    
end
    

% Now do the plots
figure; hold on;  title('Relative Errors in M - Transition Matrix'); xlabel('Num. Samples'); ylabel('Sqr. Error'); 
% Plot M errors
for i=1:HMM_x_dim
    for j=1:HMM_x_dim
       errorbar(num_points_vec, reshape(M_SQR_ERROR_MEAN(i,j,:), 1, length(num_points_vec)), reshape(M_SQR_ERROR_STD(i,j,:), 1, length(num_points_vec)));         
    end
end

figure; hold on; title('Relative Errors in MEW (divided by mew gaps) - Expectation Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error'); 
% Plot MEW errors
for i=1:HMM_x_dim
   errorbar(num_points_vec, MEW_SQR_ERROR_MEAN(i,:), MEW_SQR_ERROR_STD(i,:)); 
end

figure; hold on; title('Relative Errors in SIGMA - Standard Error Vector'); xlabel('Num. Samples'); ylabel('Sqr. Error'); 
% Plot MEW errors
for i=1:HMM_x_dim
   errorbar(num_points_vec, SIGMA_SQR_ERROR_MEAN(i,:), SIGMA_SQR_ERROR_STD(i,:)); 
end


% Plot the KL distance 
figure; hold on; title('Kullback-Leibler Error Vector'); xlabel('Num. Samples'); ylabel('KL Error'); 
% size(num_points_vec)
% size(KL_ERROR_MEAN)
% size(KL_ERROR_STD)
errorbar(num_points_vec, KL_ERROR_MEAN, KL_ERROR_STD); 


% Plot the Classification (Viterbi) error 
figure; hold on; title('Classification (Viterbi) Error Probability Vector'); xlabel('Num. Samples'); ylabel('CLASS. Error'); 
errorbar(num_points_vec, CLASS_ERROR_MEAN, CLASS_ERROR_STD);



% Plot the Bayesian Classification (forward) error 
figure; hold on; title('Classification (Forward) "Error Probability" Vector'); xlabel('Num. Samples'); ylabel('forward CLASS. Error'); 
errorbar(num_points_vec, BAYES_CLASS_ERROR_MEAN, BAYES_CLASS_ERROR_STD);