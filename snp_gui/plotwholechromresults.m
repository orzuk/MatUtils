function [DUMMY] = PlotWholeChromResults(full_path_data_file_name, full_path_model_and_output_file_name, ChromosomeToPlot, ...
                                             ExpressionPlot, ViterbiPlot, MarginalPlot, do_fold_change)


% loading the needed data 
load(full_path_data_file_name);


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


num_genes = length(HMM_genes); num_samples = length(HMM_samples); 
% First get rid of this stupid cell format
HMM_chromosome = zeros(1,num_genes);
for i=1:num_genes
    HMM_chromosome(i) = str2num(HMM_chromosome_arm{i}(find((HMM_chromosome_arm{i} <= '9')&(HMM_chromosome_arm{i} >= '0')))); % get only the numerical digits
end



plot_vec = 'rgkmc';  % Colors for plotting

% Vector denoting the start of the q-arm
first_on_q =[142604331, 95203899, 94912966, 52697919, 50705877, 62387336, 63850117, 48697951,  66551451, 41965097, 55478369, ...
        39709262, 17055261, 18814688, 18962068, 46771070, 25766888, 17462263, 34390067, 30817204, 14665357, 14608115];

epsilon = 0.000000001;


% Write Samples Nicer
for p=1:length(HMM_samples)
    HMM_samples{p}(find(HMM_samples{p} == '_')) = '-';     
end



% Now take only the samples with labels in our set (not the reference set)
% if(do_fold_change)
%     
%     HMM_take_indexes = find(HMM_ref_labels);
%     HMM_ref_indexes = find(HMM_ref_labels==0);
%     
%     % Make a row vector
%     if(size(HMM_take_indexes, 2) == 1)
%         HMM_take_indexes = HMM_take_indexes';
%         HMM_ref_indexes = HMM_ref_indexes';
%     end
%     
%     HMM_new_samples = {};
%     j=1;
%     for i=HMM_take_indexes
%         HMM_new_samples{j} = HMM_samples{i};
%         j = j + 1;       
%     end
%     HMM_samples = HMM_new_samples;
%     num_samples = length(HMM_samples);
%     

        
%     NewPatientsToPlot = [];
%     j = 1;
%     % Find the patients to plot in the new samples ... 
%     for p=PatientToPlot       
%         if(HMM_ref_labels(p) == 0)
%             print('ERRORR !!! YOU TOOK A PATIENT FROM THE REFERENCE SET !!! \n')
%         else
%             NewPatientsToPlot(j) = find(HMM_take_indexes==p);
%             j=j+1;
%         end   %if 
%         
%     end
    
%    PatientsToPlot = NewPatientsToPlot;

end

% Take all samples 
PatientsToPlot =  1:length(HMM_samples) 

% prepare legend vector 

HMM_x_dim = HMM_MODEL_Details.KMin(ChromosomeToPlot)

leg_vec = {};
if(do_fold_change)
    leg_vec{1} = 'Expression Log Ratio';
else
    leg_vec{1} = 'Expression Log';
end
end
for amp_level=HMM_x_dim:-1:2
    leg_vec{amp_level} = ['Amp. Lev. ' num2str(amp_level-1) ' Prob '];
end
leg_vec{1+HMM_x_dim} = 'Amp. Level';
leg_vec
% Done all preperations. Now do the plot : 

% Now plot the Chromosomes - for each Chromosome we plot only the average
% of amplification !!!!!
i=ChromosomeToPlot;

% First plot an overall picture of the chromosome
figure; hold on;


XMIN = 0; 
XMAX = HMM_chrom_loc_mat{i}(end)*length(HMM_samples)*1.05;
YMIN = -2; 
YMAX = HMM_x_dim;

AXIS([XMIN XMAX YMIN YMAX]);

title([ ' Whole Chrom. : ' num2str(ChromosomeToPlot)  ' (' num2str(HMM_x_dim) ' Levels)' ], 'fontsize', 8);

% Go over all patients and concatenate !!! 
HMM_whole_chrom_data = reshape(HMM_chrom_data_mat{i}, 1, size(HMM_chrom_data_mat{i}, 1) * size(HMM_chrom_data_mat{i}, 2)); % Make a vector !!
HMM_whole_chrom_Viterbi_Path = reshape(Viterbi_Path{i}, 1, size(Viterbi_Path{i}, 1) * size(Viterbi_Path{i}, 2));


HMM_whole_chrom_loc_mat = [];



GammaCumSum = {};
HMM_whole_GammaCumSum = [];
for j=PatientsToPlot
%    HMM_chrom_loc_mat{i}(end)
    % First deal with the locations 
    HMM_whole_chrom_loc_mat = [HMM_whole_chrom_loc_mat; HMM_chrom_loc_mat{i}+(j-1)*HMM_chrom_loc_mat{i}(end)];
    
   % max(HMM_whole_chrom_loc_mat)
 %   JIS = j
    
    GammaCumSum{j} = Gamma_Probs{i}(j,:,:);
    GammaCumSum{j} = cumsum(GammaCumSum{j}(:,end:-1:1,:), 2);   % Note that we do REVERSE cumsum here !!!! So we get the prob. of being at this level or higher !!! 
    HMM_whole_GammaCumSum = [ HMM_whole_GammaCumSum , reshape(GammaCumSum{j}, size(GammaCumSum{j}, 2), length(Gamma_Probs{i}(1,1,:))) ];

 %   size(reshape(GammaCumSum{j}, size(GammaCumSum{j}, 2), length(Gamma_Probs{i}(1,1,:))))
    
end
        

% Now plot the ORIGINAL expressions !!! 
if(ExpressionPlot)
    PatientToPlot=1; 
%    size(HMM_whole_chrom_loc_mat)
%    size(HMM_whole_chrom_data)
    plot( HMM_whole_chrom_loc_mat, HMM_whole_chrom_data, '*m');                                 
    
    % Give one legend to explain things 
    gap = max(max(HMM_chrom_data_mat{i}(:,:))) - min(min(HMM_chrom_data_mat{i}(:,:)));
    ylabel('Exp. Level ', 'fontsize', 8);   % Do the plot
    
end 
xlabel(['Loc. On Chrom. '  num2str(i)] , 'fontsize', 8 );


if(MarginalPlot)  % plot the marginal probabilities         
    % Compute the comulative sum
    GammaCumSum = Gamma_Probs{i}(PatientToPlot,:,:);
    GammaCumSum = cumsum(GammaCumSum(:,end:-1:1,:), 2);   % Note that we do REVERSE cumsum here !!!! So we get the prob. of being at this level or higher !!! 
    
    size(HMM_whole_GammaCumSum)
    for amp_level = 1:HMM_x_dim-1
        %                plot( HMM_chrom_loc_mat{i}, reshape(Gamma_Probs{i}(p + (f-1)*num_samples_in_page,amp_level,:),1,length(Gamma_Probs{i}(1,1,:))), ['.' plot_vec(amp_level)]); 
        plot( HMM_whole_chrom_loc_mat, HMM_whole_GammaCumSum, ['.' plot_vec(amp_level)]); 
    end
    ylabel('Amp. Level Prob. ', 'fontsize', 8);   % Do the plot
end 

if(ViterbiPlot)  % plot the most likely pathes  
    
    plot( HMM_whole_chrom_loc_mat, HMM_whole_chrom_Viterbi_Path, 'bx'); 
    ylabel('Amp. Level ', 'fontsize', 8);   % Do the plot
    
end 


if(MarginalPlot)
    legend(leg_vec);   
else
    legend(leg_vec{1:2});
end


% % Put the P-Q line
% H = line([first_on_q(i), first_on_q(i)], ...
%     [min(HMM_chrom_data_mat{i}(:,PatientToPlot)) - 0.05*gap, ...
%         max(HMM_chrom_data_mat{i}(:,PatientToPlot)) + 0.05*gap]); % set the boundary between p and q arms
% set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
% 
% H = line([first_on_q(i), first_on_q(i)], [YMIN, YMAX]); % set the boundary between p and q arms
% set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
% H = line([first_on_q(i), first_on_q(i)], [0, HMM_x_dim]); % set the boundary between p and q arms
% set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
% 

for j=PatientsToPlot
    H = line([first_on_q(i)+(j-1)*HMM_chrom_loc_mat{i}(end), first_on_q(i)+(j-1)*HMM_chrom_loc_mat{i}(end)], [YMIN, YMAX]); % set the boundary between p and q arms
    set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
    
    H = line([j*HMM_chrom_loc_mat{i}(end), j*HMM_chrom_loc_mat{i}(end)], [YMIN, YMAX]); % set the boundary between p and q arms
    set(H, 'LineWidth', 2);  set(H, 'Color', 'g');   % Color it in red 

end
    
H = line([0, XMAX], [0, 0]); % set the boundary between p and q armsset
set(H, 'LineWidth', 2);  set(H, 'Color', 'k');   % Color it in black 


figure; hist(HMM_chrom_data_mat{i}(:,PatientToPlot), 100); 
title(['Pat. : ' HMM_samples{PatientToPlot} '  Chrom. : ' num2str(ChromosomeToPlot)  'Expression Histogram' ], 'fontsize', 8);
xlabel('Log Exp. Ratio'); ylabel('Frequencey');

DUMMY = 0; 

