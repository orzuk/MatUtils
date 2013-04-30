function    DUMMY = PlotPatientOnChromResults(user_dir, SampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct);

% Various data types
MRNA_EXP = 0; SNP_CHIPS=1;
DataType = SNP_CHIPS; x_dim = HMMParamsStruct.x_dim;

hmm_out_dir = [user_dir  '/hmm_out'];
chip_type = lower(SNPChipAnnotStruct.chip);
num_samples = length(SampleNames);
genome_assembly = SNPChipAnnotStruct.genome_build;
SPECIAL_MODELS_FLAG = 1; % Must have the special SNPs model
RelevantInds = [1:6 17:22 33:38 49:54 65:70 81:86]; % This is good for the case of 3 levels

ExpressionPlot=1; MarginalPlot=1; ViterbiPlot=1; % plot everything

for(cur_sample = 1:num_samples) % Outer loop on samples ....
    sample_name = SampleNames{cur_sample};
    hmm_out_file_name = [user_dir  '/hmm_out/' sample_name '_' chip_type '_hmm_out'];
    load(hmm_out_file_name); % Load the HMM output struct

    [ALL_MAT CHROM_MATS HMMOutStruct] = LoadAndIntersectData(user_dir, sample_name, LDStruct, SNPChipAnnotStruct, HMMParamsStruct, HMMOutStruct);

    % % % % % %     % Load the needed data of the chip
    % % % % % %     load([user_dir '/' sample_name '_' chip_type '.mat']);
    % % % % % %     snp_id_str = ['snp_id_' lower(chip_type)];
    % % % % % %     eval(['snp_id_chip = ' snp_id_str ';']);
    % % % % % %     [SnpsNames I J] = intersect(SNPChipAnnotStruct.snp_ids, snp_id_chip);
    % % % % % %     SNPChipAnnotStruct.rs_ids2 = SNPChipAnnotStruct.rs_ids(I); SNPChipAnnotStruct.snp_ids2 = SNPChipAnnotStruct.snp_ids(I); % Keep also the snp ids. Do not destroy the original vec.
    % % % % % %     HMM_locations = SNPChipAnnotStruct.chr_loc_vec(I); % Currently we don't know the locations
    % % % % % %     HMM_samples = {}; HMM_samples{1} = sample_name; % Pick one of the samples
    % % % % % %     HMM_ref_labels = 1;
    % % % % % %     HMM_chromosome_arm = SNPChipAnnotStruct.chr_vec(I); % Problem ! here we've got only the chromosome and not the arm!
    % % % % % %     sample_ratio_str = ['allele_ratio_vec_' lower(chip_type)];
    % % % % % %     sample_copy_num_str = ['copy_num_vec_' lower(chip_type)];
    % % % % % %     sample_genotype_call_str = ['genotype_vec_' lower(chip_type)];
    % % % % % %     eval([sample_ratio_str '= min(' sample_ratio_str ', 9999999999);']);
    % % % % % %     eval(['HMM_data = ' sample_copy_num_str './ (' sample_ratio_str ' + 1);']);
    % % % % % %     eval(['HMM_data2 = HMM_data .* ' sample_ratio_str ';']);
    % % % % % %     HMM_data = HMM_data(J); HMM_data2 = HMM_data2(J);
    % % % % % %     TOL = 0.00000000000001;
    % % % % % %     eval(['max_error = max( (HMM_data+HMM_data2-' sample_copy_num_str ').^2)']);
    % % % % % %     SPECIAL_MODELS_FLAG=1;


    % Here's what we got :
    % HMM_chromosome_arm   - which arm (e.g. 7q) every gene lies on
    % HMM_data             - samples (expression data matrix ) of all genes
    % HMM_genes            - labels of all genes
    % HMM_locations        - locations (in nucleotides from chromosome starts of each gene on chromosome
    % HMM_samples          - labels of all samples
    % HMM_ref_labels       - labels saying if this sample is in the reference ('Normal') group or in the set we want to check ('Cancer')


    %%    load(full_path_model_and_output_file_name);

    % Here's what we got :

    % HMM_MODEL_SAVED     - The markov models
    % Viterbi_Path        - The Viterbi best pathes for each Patient and each Chromosome
    % Gamma_Probs         - The marginal probabilities for each Patient and each Chromosome
    % HMM_CHROM_KL_DIST   - The Chromosomal Distance Matrix

    % Go over all chromosomes.     % New: This part is for plotting all Chromosomes together!!!

    for use_viterbi=[0:1]
        figure; hold on; title(['All chroms SNPs Pat. : ' sample_name ' Viterbi ' num2str(use_viterbi) ]);

        ALL_MAT.chr = []; % ALL_MAT.loc = []; ALL_MAT.data = [];
        %         ALL_MAT.data_A = []; all_data_B_mat = [];
        all_viterbi_mat = []; all_viterbi_A_mat = []; all_viterbi_B_mat = []; all_viterbi_alpha_mat = []; all_viterbi_beta_mat = [];


        for i=HMMParamsStruct.ChromosomesToPlot
            do_chrom = i
            num_genes = length(CHROM_MATS.data_A{i});
            ALL_MAT.chr = [ALL_MAT.chr i+zeros(1,num_genes)];

            % Now copy the Vitebi Paths - Transfer V_VEC to its 'real' values
            do_couples=1; % flag saying whether to do couples or singletons
            joint_flag=1; % flag saying to compute marginals jointly on both chromosomes
            % %     [VitStruct.A_copylpha_genotype VitStruct.B_copyeta_genotype VitStruct.alpha_copy VitStruct.B_copyeta_copynumber ...
            % %         VitStruct.total_copy VitStruct.A_copy VitStruct.B_copy] =
            [VitStruct ProbsStruct] = GetBestMarginalPredictions(HMMOutStruct.SNPsProbsMat{i}, HMMOutStruct.ModelParams, do_couples, joint_flag);
            %            reshape(Gamma_Probs{i}(1,:,:), 256, size(Gamma_Probs{i}, 3) ), HMM_MODEL_SAVED{i}{x_dim}, do_couples);

            %%%% Use the Viterbi and not the Gamma's
% %             if(use_viterbi)
% %                 VitStruct.alpha_genotype = bitget(Viterbi_Path{i}(:,PatientToPlot),1);
% %                 VitStruct.beta_genotype = bitget(Viterbi_Path{i}(:,PatientToPlot),5);
% %                 VitStruct.alpha_copy = bitget(Viterbi_Path{i}(:,PatientToPlot),2) + ...
% %                     2*bitget(Viterbi_Path{i}(:,PatientToPlot),3) + 4*bitget(Viterbi_Path{i}(:,PatientToPlot),4);
% %                 VitStruct.beta_copy = bitget(Viterbi_Path{i}(:,PatientToPlot),6) + ...
% %                     2*bitget(Viterbi_Path{i}(:,PatientToPlot),7) + 4*bitget(Viterbi_Path{i}(:,PatientToPlot),8);
% %                 VitStruct.total_copy = VitStruct.alpha_copy+VitStruct.beta_copy;
% %                 VitStruct.A_copy = VitStruct.alpha_genotype.*VitStruct.alpha_copy + ...
% %                     VitStruct.beta_genotype.*VitStruct.beta_copy;
% %                 VitStruct.B_copy = (1-VitStruct.alpha_genotype).*VitStruct.alpha_copy + ...
% %                     (1-VitStruct.beta_genotype).*VitStruct.beta_copy;
% %             end

            all_viterbi_mat = [all_viterbi_mat VitStruct.total_copy'];
            all_viterbi_A_mat = [all_viterbi_A_mat VitStruct.A_copy'];
            all_viterbi_B_mat = [all_viterbi_B_mat VitStruct.B_copy'];
            all_viterbi_alpha_mat = [all_viterbi_alpha_mat VitStruct.alpha_copy'];
            all_viterbi_beta_mat = [all_viterbi_beta_mat VitStruct.beta_copy'];

        end
        [plot_location, chr_start_loc, end_p_location, ind] = chr_loc_into_x_axis(ALL_MAT.chr, ALL_MAT.loc, [1:22]);
        ALL_MAT.data = smooth(ALL_MAT.data, 30);
        plot(plot_location, ALL_MAT.data, 'c.');
        % Calculate the means from the data, for better displaying:
        MEAN_VEC = zeros(1,1+max(all_viterbi_mat));
        for j=1:length(MEAN_VEC)
            if(~isempty(find(all_viterbi_mat == j-1)))
                MEAN_VEC(j) = mean( ALL_MAT.data(find(all_viterbi_mat == j-1)) );
            end
        end

        load([user_dir '/' sample_name '_' chip_type '.mat']);
        snp_id_str = ['snp_id_' lower(chip_type)];
        eval(['snp_id_chip = ' snp_id_str ';']);
        [SnpsNames I J] = intersect(SNPChipAnnotStruct.snp_ids, snp_id_chip);
        sample_ratio_str = ['allele_ratio_vec_' lower(chip_type)];
        sample_copy_num_str = ['copy_num_vec_' lower(chip_type)];
        sample_genotype_call_str = ['genotype_vec_' lower(chip_type)];
        eval([sample_ratio_str '= min(' sample_ratio_str ', 9999999999);']);
        [LLL LLL_I LLL_J] = intersect(ALL_MAT.loc, SNPChipAnnotStruct.chr_loc_vec(I)); % Get the locations of the SNPs we still have

        eval(['[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(' sample_genotype_call_str '(LLL_J));']);
        %%%[no_call_ind1, AB_ind1, AA_ind1, BB_ind1] = calls_vec_into_ind(HD78_9_diag_calls_hind(LLL_J)); % Get the Genotypes
        %%%%AB_ind1 = randperm(length(LLL_J)); AB_ind1 = AB_ind1(1:13166);
        plot(plot_location, MEAN_VEC(all_viterbi_mat+1), 'k.');

        all_minmax_ratio2 = min(ALL_MAT.data_A,ALL_MAT.data_B); %%% max(ALL_MAT.data_A./ALL_MAT.data_B, ALL_MAT.data_B./ALL_MAT.data_A);
        all_minmax_ratio2 = all_minmax_ratio2(LLL_I);
        all_minmax_ratio2 = smooth(all_minmax_ratio2(AB_ind1), 30);
        all_minmax_ratio = min(ALL_MAT.data_A,ALL_MAT.data_B)./max(ALL_MAT.data_A,ALL_MAT.data_B);
        all_minmax_ratio = all_minmax_ratio(LLL_I);
        all_minmax_ratio = smooth(all_minmax_ratio, 30);
        all_min = min(ALL_MAT.data_A,ALL_MAT.data_B);
        all_max = max(ALL_MAT.data_A,ALL_MAT.data_B);
        all_min = smooth(all_min, 30);
        all_max = smooth(all_max, 30);
        ALL_MAT.data_A = smooth(ALL_MAT.data_A, 30);
        ALL_MAT.data_B = smooth(ALL_MAT.data_B, 30);


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

    num_genes = length(HMM_genes);
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
    for p=1:num_samples
        sample_name(find(sample_name == '_')) = '-';
    end

    leg_vec = {};
    if(do_fold_change)
        leg_vec{1} = 'Expression Log Ratio';
    else
        leg_vec{1} = 'Expression Log';
    end
    for amp_level=x_dim:-1:2
        leg_vec{amp_level} = ['Amp. Lev. ' num2str(amp_level-1) ' Prob '];
    end
    leg_vec{1+x_dim} = 'Amp. Level';
    leg_vec
    % Done all preperations. Now do the plot :

    % Now plot the Chromosomes - for each Chromosome we plot only the average
    % of amplification !!!!!
    i=ChromosomeToPlot;

    % First plot an overall picture of the chromosome
    figure; hold on;


    XMIN = 0; XMAX = HMM_chrom_loc_mat{i}(end)*1.05; YMIN = -0.2; YMAX = 2*x_dim-1;
    AXIS([XMIN XMAX YMIN YMAX]);

    title(['Pat. : ' sample_name '  Chrom. : ' num2str(ChromosomeToPlot)  ' (' num2str(x_dim) ...
        ' Levels)  Mean Exp. ' num2str(mean(CHROMS_MAT.data{i})) ' Std. Exp. '  num2str(std(CHROMS_MAT.data{i}))], 'fontsize', 8);
    % Now plot the ORIGINAL expressions !!!
    if(ExpressionPlot)
        do_clip=0;
        if(do_clip)
            EPS = 0.1;
            CLIPPED_VEC = ClipModelObservedData( HMM_MODEL_SAVED{ChromosomeToPlot}{2}, CHROMS_MAT.data{i}) + EPS;
        else
            CLIPPED_VEC = CHROMS_MAT.data{i}(:,PatientToPlot);
            if(SPECIAL_MODELS_FLAG==1)
                CLIPPED_VEC2 = CHROMS_MAT.data2{i}(:,PatientToPlot);
            end
        end  % do clip

        plot( HMM_chrom_loc_mat{i}, CLIPPED_VEC, '*m');
        if(SPECIAL_MODELS_FLAG==1)
            plot( HMM_chrom_loc_mat{i}, CLIPPED_VEC2, '*g');
        end

        % Give one legend to explain things
        gap = max(CHROMS_MAT.data{i}(:,PatientToPlot)) - min(CHROMS_MAT.data{i}(:,PatientToPlot));
        ylabel('Exp. Level ', 'fontsize', 8);   % Do the plot

    end
    xlabel(['Loc. On Chrom. '  num2str(i)] , 'fontsize', 8 );

    if(MarginalPlot)  % plot the marginal probabilities
        % Compute the comulative sum
        GammaCumSum = Gamma_Probs{i}(PatientToPlot,:,:);
        GammaCumSum = cumsum(GammaCumSum(:,end:-1:1,:), 2);   % Note that we do REVERSE cumsum here !!!! So we get the prob. of being at this level or higher !!!

        for amp_level = 1:x_dim-1
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
        [min(CHROMS_MAT.data{i}(:,PatientToPlot)) - 0.05*gap, ...
        max(CHROMS_MAT.data{i}(:,PatientToPlot)) + 0.05*gap]); % set the boundary between p and q arms
    set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red

    H = line([first_on_q(i), first_on_q(i)], [YMIN, YMAX]); % set the boundary between p and q arms
    set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red
    H = line([first_on_q(i), first_on_q(i)], [0, x_dim]); % set the boundary between p and q arms
    set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red

    H = line([0, XMAX], [0, 0]); % set the boundary between p and q armsset
    set(H, 'LineWidth', 2);  set(H, 'Color', 'k');   % Color it in red


    figure; hist(CHROMS_MAT.data{i}(:,PatientToPlot), 100);
    title(['Pat. : ' sample_name '  Chrom. : ' num2str(ChromosomeToPlot)  ' Expression Histogram' ], 'fontsize', 8);
    xlabel('Log Exp. Ratio'); ylabel('Frequencey');

end % loop on patients

DUMMY = 0;



