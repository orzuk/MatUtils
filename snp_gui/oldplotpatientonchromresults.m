function [DUMMY] = PlotPatientOnChromResults(full_path_data_file_name, full_path_model_and_output_file_name, DataType, ...
    ChromosomeToPlot, PatientToPlot, ...
    ExpressionPlot, ViterbiPlot, MarginalPlot, do_fold_change, do_clip)

% Various data types
MRNA_EXP = 0; SNP_CHIPS=1;

% loading the needed data 
load(full_path_data_file_name);

% In the new version, we  need to do some conversions of the data ...
if(DataType == SNP_CHIPS) % Work for now only on the Hind SNP chip
    load('LD_all_chroms.mat')
    load('Hind_allele_info_table.mat'); % Load the Hind details

    [SnpsNames I J] = intersect(hind_table(2:end,1), snp_id_hind);

    HMM_genes = snp_id_hind;
    HMM_locations = str2num(char(hind_table(I,3))); % Currently we don't know the locations
    HMM_samples = {}; HMM_samples{1} = 'HD78_9_diag'; % Pick one of the samples
    sample_ratio_str = [HMM_samples{1} '_ratio_vec_hind'];
    sample_copy_num_str = [HMM_samples{1} '_raw_copy_num_hind'];
    sample_genotype_call_str = [HMM_samples{1} '_calls_hind'];
    HMM_ref_labels = 1;
    HMM_chromosome_arm = hind_table(I,2); % Problem ! here we've got only the chromosome and not the arm!
    eval([sample_ratio_str '= min(' sample_ratio_str ', 9999999999);']);
    eval(['HMM_data = ' sample_copy_num_str './ (' sample_ratio_str ' + 1);']);
    eval(['HMM_data2 = HMM_data .* ' sample_ratio_str ';']);
% % % %     HD78_9_diag_ratio_vec_hind = min(HD78_9_diag_ratio_vec_hind, 9999999999); % Avoid infinities
% % % %     HMM_data = HD78_9_diag_raw_copy_num_hind ./ (HD78_9_diag_ratio_vec_hind+1);
% % % %     HMM_data2 = HMM_data .* HD78_9_diag_ratio_vec_hind;
    HMM_data = HMM_data(J); HMM_data2 = HMM_data2(J);
    SPECIAL_MODELS_FLAG=1;

end

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

% Go over all chromosomes.     % New: This part is for plotting all Chromosomes together!!!
  
for use_viterbi=[0:1]
figure; hold on; title(['All chroms SNPs Pat. : ' HMM_samples{1} ' Viterbi ' num2str(use_viterbi) ]);

all_chrom_mat = []; all_chrom_loc_mat = []; all_data_mat = []; all_viterbi_mat = [];
all_data_A_mat = []; all_data_B_mat = [];
all_viterbi_A_mat = []; all_viterbi_B_mat = []; all_viterbi_alpha_mat = []; all_viterbi_beta_mat = [];

HMM_x_dim = HMM_MODEL_Details.KMin(ChromosomeToPlot)
for i=1:22 %ChromosomeToPlot
    do_chrom = i
    num_genes = length(HMM_chrom_loc_mat{i});
    all_chrom_mat = [all_chrom_mat i+zeros(1,num_genes)];
    all_chrom_loc_mat = [all_chrom_loc_mat HMM_chrom_loc_mat{i}'];
    all_data_mat = [all_data_mat HMM_chrom_data_mat{i}'+HMM_chrom_data_mat2{i}'];
    all_data_A_mat = [all_data_A_mat HMM_chrom_data_mat{i}'];
    all_data_B_mat = [all_data_B_mat HMM_chrom_data_mat2{i}'];
    size(Gamma_Probs{i}, 3)
    % Now copy the Vitebi Paths - Transfer V_VEC to its 'real' values
    do_couples=1; % flag saying whether to do couples or singletons
% %     [VitStruct.A_copylpha_genotype VitStruct.B_copyeta_genotype VitStruct.alpha_copy VitStruct.B_copyeta_copynumber ...
% %         VitStruct.total_copy VitStruct.A_copy VitStruct.B_copy] =     
     [VitStruct ProbsStruct] = GetBestMarginalPredictions(reshape(Gamma_Probs{i}(1,:,:), 256, size(Gamma_Probs{i}, 3) ), HMM_MODEL_SAVED{i}{HMM_x_dim}, do_couples);

    %%%% Use the Viterbi and not the Gamma's
    if(use_viterbi)
        VitStruct.alpha_genotype = bitget(Viterbi_Path{i}(:,PatientToPlot),1);
        VitStruct.beta_genotype = bitget(Viterbi_Path{i}(:,PatientToPlot),5);
        VitStruct.alpha_copy = bitget(Viterbi_Path{i}(:,PatientToPlot),2) + ...
            2*bitget(Viterbi_Path{i}(:,PatientToPlot),3) + 4*bitget(Viterbi_Path{i}(:,PatientToPlot),4);
        VitStruct.beta_copy = bitget(Viterbi_Path{i}(:,PatientToPlot),6) + ...
            2*bitget(Viterbi_Path{i}(:,PatientToPlot),7) + 4*bitget(Viterbi_Path{i}(:,PatientToPlot),8);
        VitStruct.total_copy = VitStruct.alpha_copy+VitStruct.beta_copy;
        VitStruct.A_copy = VitStruct.alpha_genotype.*VitStruct.alpha_copy + ...
            VitStruct.beta_genotype.*VitStruct.beta_copy;
        VitStruct.B_copy = (1-VitStruct.alpha_genotype).*VitStruct.alpha_copy + ...
            (1-VitStruct.beta_genotype).*VitStruct.beta_copy;
    end

    all_viterbi_mat = [all_viterbi_mat VitStruct.total_copy'];
    all_viterbi_A_mat = [all_viterbi_A_mat VitStruct.A_copy'];
    all_viterbi_B_mat = [all_viterbi_B_mat VitStruct.B_copy'];
    all_viterbi_alpha_mat = [all_viterbi_alpha_mat VitStruct.alpha_copy'];
    all_viterbi_beta_mat = [all_viterbi_beta_mat VitStruct.beta_copy'];

end
[plot_location, chr_start_loc, end_p_location, ind] = chr_loc_into_x_axis(all_chrom_mat, all_chrom_loc_mat, [1:22]);
all_data_mat = smooth(all_data_mat, 30);
plot(plot_location, all_data_mat, 'c.');
% Calculate the means from the data, for better displaying:

MEAN_VEC = zeros(1,1+max(all_viterbi_mat));
for j=1:length(MEAN_VEC)
    if(~isempty(find(all_viterbi_mat == j-1)))
       MEAN_VEC(j) = mean( all_data_mat(find(all_viterbi_mat == j-1)) );
    end
end


[LLL LLL_I LLL_J] = intersect(all_chrom_loc_mat, HMM_locations); % Get the locations of the SNPs we still have

eval(['[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(' sample_genotype_call_str '(LLL_J));']);
%%%[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(HD78_9_diag_calls_hind(LLL_J)); % Get the Genotypes 
%%%%AB_ind1 = randperm(length(LLL_J)); AB_ind1 = AB_ind1(1:13166);
plot(plot_location, MEAN_VEC(all_viterbi_mat+1), 'k.');

all_minmax_ratio2 = min(all_data_A_mat,all_data_B_mat); %%% max(all_data_A_mat./all_data_B_mat, all_data_B_mat./all_data_A_mat);
all_minmax_ratio2 = all_minmax_ratio2(LLL_I);
all_minmax_ratio2 = smooth(all_minmax_ratio2(AB_ind1), 30);
all_minmax_ratio = min(all_data_A_mat,all_data_B_mat)./max(all_data_A_mat,all_data_B_mat);
all_minmax_ratio = all_minmax_ratio(LLL_I);
all_minmax_ratio = smooth(all_minmax_ratio, 30);
all_min = min(all_data_A_mat,all_data_B_mat); 
all_max = max(all_data_A_mat,all_data_B_mat); 
all_min = smooth(all_min, 30);
all_max = smooth(all_max, 30);
all_data_A_mat = smooth(all_data_A_mat, 30);
all_data_B_mat = smooth(all_data_B_mat, 30);



% Plot the 'data' of A and B seperately
plot(plot_location, all_min-4, 'y.');
plot(plot_location, all_max-8, 'c.');

%%plot(plot_location(LLL_I(AB_ind1)), (all_minmax_ratio2)-4, 'm.');
%plot(plot_location(LLL_I(AB_ind1)), 1./(all_minmax_ratio(AB_ind1))-8, 'c.');

%%plot(plot_location, log2(1./all_minmax_ratio)-8, 'c.');

% Plot the four 'partial' genomes
%plot(plot_location, all_viterbi_A_mat-12, 'r');
%plot(plot_location, all_viterbi_B_mat-16, 'g');
plot(plot_location, all_viterbi_alpha_mat-12, 'r');
plot(plot_location, all_viterbi_beta_mat-16, 'b');

% Count how many 'switches' we have:
A_switches = sum(all_viterbi_A_mat(1:end-1) ~= all_viterbi_A_mat(2:end))
B_switches = sum(all_viterbi_B_mat(1:end-1) ~= all_viterbi_B_mat(2:end))
alpha_switches = sum(all_viterbi_alpha_mat(1:end-1) ~= all_viterbi_alpha_mat(2:end))
beta_switches = sum(all_viterbi_beta_mat(1:end-1) ~= all_viterbi_beta_mat(2:end))

% Plot the chromosomes start and end lines
chr_vec = 1:22;
ax = axis;
for j=1:length(chr_start_loc)
    plot([chr_start_loc(j) chr_start_loc(j)], ax(3:4),'k');
end
plot(ax(1:2),[0 0],'k');
set(gca, 'xtick',chr_start_loc,'xticklabel',chr_vec, 'FontSize', 7);
for j=1:length(end_p_location)
    plot([end_p_location(j) end_p_location(j)], ax(3:4),'-.k');
end

YLims = get(gca, 'YLim');
set(gca, 'ytick',[YLims(1):YLims(2)],'yticklabel',[0:3 0:3 0:3 0:3 0:4]);
legend('raw copy', 'mean', 'min copy', 'max copy', '\alpha', '\beta');
end % For on Viterbi 

DUMMY = 0; 
return; % drop all that is below

num_genes = length(HMM_genes); num_samples = length(HMM_samples); 
% First get rid of this stupid cell format
% Get rid of this stupid cell format, and get data as numbers
if(DataType == SNP_CHIPS)
    HMM_chromosome = (char(HMM_chromosome_arm)); % Here we already have the chromosomes only
    XXX = find(HMM_chromosome == 'X');
    HMM_chromosome(XXX,1) = '2'; % Encode the X chromosome
    HMM_chromosome(XXX,2) = '3'; % Encode the X chromosome
    III = intersect( find(HMM_chromosome(:,1) == ' '), find(HMM_chromosome(:,2) == ' ') );
    HMM_chromosome(III,1) = '0'; % Encode the chromosomes not known at all
    HMM_chromosome = str2num(HMM_chromosome);
else
    HMM_chromosome = zeros(1,num_genes);
    for i=1:num_genes
        HMM_chromosome(i) = str2num(HMM_chromosome_arm{i}(find((HMM_chromosome_arm{i} <= '9')&(HMM_chromosome_arm{i} >= '0')))); % get only the numerical digits
    end
end


plot_vec = 'rgkmc';  % Colors for plotting

% Vector denoting the start of the q-arm
first_on_q =[142604331, 95203899, 94912966, 52697919, 50705877, 62387336, 63850117, 48697951,  66551451, 41965097, 55478369, ...
        39709262, 17055261, 18814688, 18962068, 46771070, 25766888, 17462263, 34390067, 30817204, 14665357, 14608115];

epsilon = 0.000000001;

num_chroms = 22;


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
%     
% end

% prepare legend vector 


leg_vec = {};
if(do_fold_change)
    leg_vec{1} = 'Expression Log Ratio';
else
    leg_vec{1} = 'Expression Log';
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
XMAX = HMM_chrom_loc_mat{i}(end)*1.05;
YMIN = -0.2; 
YMAX = 2*HMM_x_dim-1;

AXIS([XMIN XMAX YMIN YMAX]);

title(['Pat. : ' HMM_samples{PatientToPlot} '  Chrom. : ' num2str(ChromosomeToPlot)  ' (' num2str(HMM_x_dim) ...
        ' Levels)  Mean Exp. ' num2str(mean(HMM_chrom_data_mat{i}(:,PatientToPlot))) ' Std. Exp. '  num2str(std(HMM_chrom_data_mat{i}(:,PatientToPlot)))], 'fontsize', 8);
% Now plot the ORIGINAL expressions !!! 
if(ExpressionPlot)
    
    if(do_clip)
       EPS = 0.1;
       CLIPPED_VEC = ClipModelObservedData( HMM_MODEL_SAVED{ChromosomeToPlot}{2}, HMM_chrom_data_mat{i}(:,PatientToPlot)) + EPS; 
    else
       CLIPPED_VEC = HMM_chrom_data_mat{i}(:,PatientToPlot);
       if(SPECIAL_MODELS_FLAG==1)     
           CLIPPED_VEC2 = HMM_chrom_data_mat2{i}(:,PatientToPlot);
       end
    end  % do clip
    
    plot( HMM_chrom_loc_mat{i}, CLIPPED_VEC, '*m');      
    if(SPECIAL_MODELS_FLAG==1) 
       plot( HMM_chrom_loc_mat{i}, CLIPPED_VEC2, '*g');                                                                
    end
       
    % Give one legend to explain things 
    gap = max(HMM_chrom_data_mat{i}(:,PatientToPlot)) - min(HMM_chrom_data_mat{i}(:,PatientToPlot));
    ylabel('Exp. Level ', 'fontsize', 8);   % Do the plot
    
end 
xlabel(['Loc. On Chrom. '  num2str(i)] , 'fontsize', 8 );



if(MarginalPlot)  % plot the marginal probabilities         
    % Compute the comulative sum
    GammaCumSum = Gamma_Probs{i}(PatientToPlot,:,:);
    GammaCumSum = cumsum(GammaCumSum(:,end:-1:1,:), 2);   % Note that we do REVERSE cumsum here !!!! So we get the prob. of being at this level or higher !!! 
    
    for amp_level = 1:HMM_x_dim-1
        %                plot( HMM_chrom_loc_mat{i}, reshape(Gamma_Probs{i}(p + (f-1)*num_samples_in_page,amp_level,:),1,length(Gamma_Probs{i}(1,1,:))), ['.' plot_vec(amp_level)]); 
        plot( HMM_chrom_loc_mat{i}, reshape(GammaCumSum(1,amp_level,:),1,length(Gamma_Probs{i}(1,1,:))), ['.' plot_vec(amp_level)]); 
    end
    ylabel('Amp. Level Prob. ', 'fontsize', 8);   % Do the plot
end 

if(ViterbiPlot)  % plot the most likely pathes  
    if(DataType == SNP_CHIPS) % Work for now only on the Hind SNP chip
        % Transfer X_VEC to its 'real' values
        VitStruct.alpha_genotype = bitget(Viterbi_Path{i}(:,PatientToPlot),1);
        VitStruct.beta_genotype = bitget(Viterbi_Path{i}(:,PatientToPlot),5);
        VitStruct.alpha_copy = bitget(Viterbi_Path{i}(:,PatientToPlot),2) + ...
            2*bitget(Viterbi_Path{i}(:,PatientToPlot),3) + 4*bitget(Viterbi_Path{i}(:,PatientToPlot),4);
        VitStruct.beta_copy = bitget(Viterbi_Path{i}(:,PatientToPlot),6) + ...
            2*bitget(Viterbi_Path{i}(:,PatientToPlot),7) + 4*bitget(Viterbi_Path{i}(:,PatientToPlot),8);
        VitStruct.total_copy = VitStruct.alpha_copy+VitStruct.beta_copy;
        VitStruct.A_copy = VitStruct.alpha_genotype.*VitStruct.alpha_copy + ...
            VitStruct.beta_genotype.*VitStruct.beta_copy;
        VitStruct.B_copy = (1-VitStruct.alpha_genotype).*VitStruct.alpha_copy + ...
            (1-VitStruct.beta_genotype).*VitStruct.beta_copy;
        
        plot(HMM_chrom_loc_mat{i}, VitStruct.A_copy+0.05, 'b.'); 
        plot(HMM_chrom_loc_mat{i}, VitStruct.B_copy-0.05, 'r.');
        plot(HMM_chrom_loc_mat{i}, VitStruct.total_copy, 'c.');
%         plot(HMM_chrom_loc_mat{i}, VitStruct.alpha_copy+0.02, 'k.');
%         plot(HMM_chrom_loc_mat{i}, VitStruct.beta_copy-0.02, 'g.');
    else    
        plot( HMM_chrom_loc_mat{i}, Viterbi_Path{i}(PatientToPlot,:), 'bx'); 
    end
    ylabel('Amp. Level ', 'fontsize', 8);   % Do the plot    
end


if(MarginalPlot)
    legend(leg_vec);   
else
    legend(leg_vec{1:2});
end
% Put the P-Q line
H = line([first_on_q(i), first_on_q(i)], ...
    [min(HMM_chrom_data_mat{i}(:,PatientToPlot)) - 0.05*gap, ...
        max(HMM_chrom_data_mat{i}(:,PatientToPlot)) + 0.05*gap]); % set the boundary between p and q arms
set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 

H = line([first_on_q(i), first_on_q(i)], [YMIN, YMAX]); % set the boundary between p and q arms
set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 
H = line([first_on_q(i), first_on_q(i)], [0, HMM_x_dim]); % set the boundary between p and q arms
set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red 

H = line([0, XMAX], [0, 0]); % set the boundary between p and q armsset
set(H, 'LineWidth', 2);  set(H, 'Color', 'k');   % Color it in red 


figure; hist(HMM_chrom_data_mat{i}(:,PatientToPlot), 100); 
title(['Pat. : ' HMM_samples{PatientToPlot} '  Chrom. : ' num2str(ChromosomeToPlot)  ' Expression Histogram' ], 'fontsize', 8);
xlabel('Log Exp. Ratio'); ylabel('Frequencey');

DUMMY = 0; 