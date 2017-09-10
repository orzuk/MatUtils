% A function for learning the models from the data. Here is what the
% program can do, according to the WhatToDo parameter which have one bit
% for every option :
% lsb - If we want to learn model parameters
% 2nd lsb - If we want to output Viterbi pathes for all patients
% 3rd lsb - If we want to output Gamma probs for all patients
% 4th lsb - If we want to output the Chromosomal Distance Matrix
% 5th lsb - If to do center&normalize or to do fold-change with respect to normal tissues

function [DUMMY] = PlotHMMChromResults(full_path_data_file_name, full_path_model_and_output_file_name, HMM_x_dim, ExpressionPlot, ViterbiPlot, MarginalPlot, ...
    ChromosomesToPlot, PatientsToPlot, ChromosomesDistancePlot, ...
    num_samples_in_row, num_samples_in_column, do_fold_change)

% load data from some chromosome
path(path,'E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom');  % weizmann
path(path,'C:\Weizmann\HMM_ido_2004_02_16\hmm_chromosome\hmm_chrom\hmm_chrom');  % home
path(path,'C:\Weizmann\HMM_ido_2004_02_16\hmm_chromosome\hmm_chrom');  % home


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
    %    HMM_chromosome(i) = str2num(HMM_chromosome_arm{i}(1:end-1));
    HMM_chromosome(i) = str2num(HMM_chromosome_arm{i}(find((HMM_chromosome_arm{i} <= '9')&(HMM_chromosome_arm{i} >= '0')))); % get only the numerical digits
end

plot_vec = 'byrgmck';  % Colors for plotting

% Vector denoting the start of the q-arm
first_on_q =[142604331, 95203899, 94912966, 52697919, 50705877, 62387336, 63850117, 48697951,  66551451, 41965097, 55478369, ...
    39709262, 17055261, 18814688, 18962068, 46771070, 25766888, 17462263, 34390067, 30817204, 14665357, 14608115];

epsilon = 0.000000001;
num_chroms = 22;

%  plot the Chromosomal Distance Matrix
if(ChromosomesDistancePlot)
    figure; hold on; imagesc(HMM_CHROM_KL_DIST); colorbar;
    title('Chromosomes Distance Matrix Between Models'); xlabel('chrom. num'); ylabel('chrom. num.');
end

% Write Samples Nicer
for p=1:length(HMM_samples)
    HMM_samples{p}(find(HMM_samples{p} == '_')) = '-';
end

% Now take only the samples with labels in our set (not the reference set)
if(do_fold_change)
    HMM_take_indexes = find(HMM_ref_labels);
    HMM_ref_indexes = find(HMM_ref_labels==0);

    % Make a row vector
    if(size(HMM_take_indexes, 2) == 1)
        HMM_take_indexes = HMM_take_indexes';
        HMM_ref_indexes = HMM_ref_indexes';
    end
    HMM_new_samples = {};
    j=1;
    for i=HMM_take_indexes
        HMM_new_samples{j} = HMM_samples{i};
        j = j + 1;
    end
    HMM_samples = HMM_new_samples;
    num_samples = length(HMM_samples);

    NewPatientsToPlot = [];
    j = 1;
    % Find the patients to plot in the new samples ...
    for p=PatientsToPlot
        if(HMM_ref_labels(p) == 0)
            print('ERRORR !!! YOU TOOK A PATIENT FROM THE REFERENCE SET !!! \n')
        else
            NewPatientsToPlot(j) = find(HMM_take_indexes==p);
            j=j+1;
        end   %if
    end
    PatientsToPlot = NewPatientsToPlot;
end

% prepare legend vector
leg_vec = {};
for amp_level=HMM_x_dim-1:-1:1
    leg_vec{amp_level} = ['amp. lev. ' num2str(amp_level)];
end


% Now plot the patients - each patient has two pages !!!!!!
for p=PatientsToPlot

    % Plot the Viterbi Pathes
    if(ViterbiPlot)  % plot the most likely pathes
        for ploty=0:1         % plot first Chromosomes 1-11 and then 12-22
            figure; hold on; title(['Patient : ' num2str(p) ' is : ' HMM_samples{p} ' Page : ' num2str(ploty+1)], 'fontsize', 8); subplot(3,4,1); hold on; title(['Patient : ' num2str(p) ' is : ' HMM_samples{p} ' Page : ' num2str(ploty+1)], 'fontsize', 8);
            for i=1:11
                subplot(3,4,i); hold on;
                p_IS = p
                INNN = i+11*ploty
                locsize = size(HMM_chrom_loc_mat{i+11*ploty})
                vitsize = size(Viterbi_Path{i+11*ploty})
                plot( HMM_chrom_loc_mat{i+11*ploty}, Viterbi_Path{i+11*ploty}(p,:), '.');  ylabel('Amp. Level', 'fontsize', 8);   % Do the plot

                H = line([first_on_q(i+11*ploty), first_on_q(i+11*ploty)], [0, HMM_x_dim-1]); % set the boundary between p and q arms
                set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red
                if(i == 10)
                    xlabel(['Loc. Chrom. '], 'fontsize', 8 );
                end
                title(['Chrom. ' num2str(i+11*ploty) ], 'fontsize', 8);  %%% 'Chromosome ' num2str(i) ' with  ' num2str(HMM_x_dim) ' levels ' num2str(y_dim) ' mixtures'] );
            end
            subplot(3,4,12); hold on;
            text(0.2,0.8,['Patient # : ' num2str(p)]); text(0.2,0.2, [' Page : ' num2str(ploty+1)], 'fontsize', 8); text(0.2,0.5, ['Name : ' HMM_samples{p}], 'fontsize', 8);
        end
    end


    % Plot the Gamma Probes
    if(MarginalPlot)
        for ploty=0:1         % plot first Chromosomes 1-11 and then 12-22
            figure; hold on; title(['Patient : ' num2str(p) ' is : ' HMM_samples{p} ' Page : ' num2str(ploty+1)], 'fontsize', 8); subplot(3,4,1); hold on; title(['Patient : ' num2str(p) ' is : ' HMM_samples{p} ' Page : ' num2str(ploty+1)], 'fontsize', 8);
            for i=1:11



                subplot(3,4,i); hold on;

                if(HMM_MODEL_Details.KMin(i+11*ploty) > 1)   % This means that we have more than one level !!!!!

                    % Compute the comulative sum
                    GammaCumSum = Gamma_Probs{i+11*ploty}(p,:,:);
                    GammaCumSum = cumsum(GammaCumSum(:,end:-1:1,:), 2);   % Note that we do REVERSE cumsum here !!!! So we get the prob. of being at this level or higher !!!



                    %     Gamma_Probs_Size = size(Gamma_Probs{i+11*ploty})
                    %     GammaCumSum_Size = size(GammaCumSum)

                    % Now plot all the probs
                    for amp_level = 1:HMM_MODEL_Details.KMin(i+11*ploty)-1  %% HMM_x_dim-1
                        %           HMM_chrom_loc_mat{i+11*ploty}
                        %           GammaCumSum(1,amp_level,:)
                        %           index_is = i+11*ploty
                        %           dim_is = HMM_MODEL_Details.KMin(i+11*ploty)
                        %           p
                        %           size(Gamma_Probs{i+11*ploty})
                        %           Gamma_Probs{i+11*ploty}(p,2,:)
                        %           plot_vec(amp_level)

                        %%%                    plot( HMM_chrom_loc_mat{i+11*ploty}, reshape(Gamma_Probs{i+11*ploty}(p,2,:), 1, length(Gamma_Probs{i+11*ploty}(p,2,:))), '.');
                        plot( HMM_chrom_loc_mat{i+11*ploty}, reshape(GammaCumSum(1,amp_level,:), 1, length(Gamma_Probs{i+11*ploty}(p,2,:))), ['.' plot_vec(amp_level)]);


                        %%%%                  plot( HMM_chrom_loc_mat{i}, reshape(GammaCumSum(1,amp_level,:),1,length(Gamma_Probs{i}(1,1,:))), ['.' plot_vec(amp_level)]);
                    end
                    ylabel('Amp. Level Prob.', 'fontsize', 8);   % Do the plot

                    % Give one legend to explain things
                    if(i==1)
                        legend(leg_vec);
                    end


                    H = line([first_on_q(i+11*ploty), first_on_q(i+11*ploty)], [0, 1]); % set the boundary between p and q arms
                    set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red
                    if(i == 10)
                        xlabel(['Loc. Chrom. '], 'fontsize', 8 );
                    end
                    title(['Chrom. ' num2str(i+11*ploty) ], 'fontsize', 8);  %%% 'Chromosome ' num2str(i) ' with  ' num2str(HMM_x_dim) ' levels ' num2str(y_dim) ' mixtures'] );
                end
                subplot(3,4,12); hold on;
                text(0.2,0.8,['Patient # : ' num2str(p)]); text(0.2,0.2, [' Page : ' num2str(ploty+1)], 'fontsize', 8); text(0.2,0.5, ['Name : ' HMM_samples{p}], 'fontsize', 8);
            end % if more than one level
        end   % loop on i
    end

end





% Now plot the Chromosomes - for each Chromosome we plot only the average
% of amplification !!!!!
for i=ChromosomesToPlot

    %  size(Gamma_Probs{i})

    % First plot an overall picture of the chromosome
    figure; hold on;



    Gamma_mean{i} = zeros(1,length(Gamma_Probs{i}(1,1,:)));
    Gamma_std{i} = zeros(1,length(Gamma_Probs{i}(1,1,:)));
    for p=1:num_samples
        for amp_level = 1:HMM_x_dim
            %    Gamma_mean{i} = Gamma_mean{i} + Viterbi_Path{i}(p,:);
            Gamma_mean{i} = Gamma_mean{i} + (amp_level-1)*reshape(Gamma_Probs{i}(p,amp_level,:),1,length(Gamma_Probs{i}(1,1,:)));
            Gamma_std{i} = Gamma_std{i} + (amp_level-1)^2*(reshape(Gamma_Probs{i}(p,amp_level,:),1,length(Gamma_Probs{i}(1,1,:))));
        end
    end
    Gamma_mean{i} = Gamma_mean{i} ./ num_samples;
    Gamma_std{i} = Gamma_std{i} ./ num_samples;
    Gamma_std{i} = Gamma_std{i} - (Gamma_mean{i}).^2;
    title(['Chromosome : ' num2str(i) ' Average Amplification Rate (with srd.) : '], 'fontsize', 8);
    errorbar( HMM_chrom_loc_mat{i}, Gamma_mean{i}, Gamma_std{i}, '.');  ylabel('Amp. Level Average', 'fontsize', 8);   % Do the plot
    H = line([first_on_q(i), first_on_q(i)], [0, HMM_x_dim-1]); % set the boundary between p and q arms
    set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red






    % Now we plot For ALL the patients this Chromosome !!
    num_samples_in_page = num_samples_in_row * num_samples_in_column;
    num_figures = ceil(num_samples/num_samples_in_page)





    % Now plot the ORIGINAL expressions !!!
    if(ExpressionPlot)
        for f=1:num_figures
            figure;
            for p=1:min(num_samples_in_page, num_samples-(f-1)*num_samples_in_page)

                subplot(num_samples_in_column, num_samples_in_row, p); hold on;

                plot( HMM_chrom_loc_mat{i}, HMM_chrom_data_mat{i}(:,p + (f-1)*num_samples_in_page), '.');

                % Give one legend to explain things
                %  if(p==1)
                %      legend(leg_vec);
                %  end
                gap = max(HMM_chrom_data_mat{i}(:,p + (f-1)*num_samples_in_page)) - min(HMM_chrom_data_mat{i}(:,p + (f-1)*num_samples_in_page));
                ylabel('Exp. Level ', 'fontsize', 8);   % Do the plot
                H = line([first_on_q(i), first_on_q(i)], ...
                    [min(HMM_chrom_data_mat{i}(:,p + (f-1)*num_samples_in_page)) - 0.05*gap, ...
                    max(HMM_chrom_data_mat{i}(:,p + (f-1)*num_samples_in_page)) + 0.05*gap]); % set the boundary between p and q arms
                set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red
                title(['Pat. : ' HMM_samples{p}], 'fontsize', 8);

                if(p == min(num_samples_in_page, num_samples-(f-1)*num_samples_in_page))
                    xlabel(['Loc. Chrom. '  num2str(i)] , 'fontsize', 8 );
                end

            end
        end
    end




    if(ViterbiPlot)  % plot the most likely pathes
        for f=1:num_figures
            figure;
            for p=1:min(num_samples_in_page, num_samples-(f-1)*num_samples_in_page)
                %%%cur_fig = mod(p, num_samples_in_page)+1;
                subplot(num_samples_in_column, num_samples_in_row, p); hold on;

                plot( HMM_chrom_loc_mat{i}, Viterbi_Path{i}(p + (f-1)*num_samples_in_page,:), '.');


                % Give one legend to explain things
                %  if(p==1)
                %      legend(leg_vec);
                %  end
                ylabel('Amp. Level ', 'fontsize', 8);   % Do the plot
                H = line([first_on_q(i), first_on_q(i)], [0, HMM_x_dim]); % set the boundary between p and q arms
                set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red
                title(['Pat. : ' HMM_samples{p}], 'fontsize', 8);

                if(p == min(num_samples_in_page, num_samples-(f-1)*num_samples_in_page))
                    xlabel(['Loc. Chrom. '  num2str(i)] , 'fontsize', 8 );
                end

            end
        end
    end



    if(MarginalPlot)  % plot the marginal probabilities
        for f=1:num_figures
            figure;
            for p=1:min(num_samples_in_page, num_samples-(f-1)*num_samples_in_page)
                %%%cur_fig = mod(p, num_samples_in_page)+1;
                subplot(num_samples_in_column, num_samples_in_row, p); hold on;


                % Compute the comulative sum
                GammaCumSum = Gamma_Probs{i}(p + (f-1)*num_samples_in_page,:,:);
                GammaCumSum = cumsum(GammaCumSum(:,end:-1:1,:), 2);   % Note that we do REVERSE cumsum here !!!! So we get the prob. of being at this level or higher !!!

                for amp_level = 1:HMM_x_dim-1
                    %                plot( HMM_chrom_loc_mat{i}, reshape(Gamma_Probs{i}(p + (f-1)*num_samples_in_page,amp_level,:),1,length(Gamma_Probs{i}(1,1,:))), ['.' plot_vec(amp_level)]);
                    plot( HMM_chrom_loc_mat{i}, reshape(GammaCumSum(1,amp_level,:),1,length(Gamma_Probs{i}(1,1,:))), ['.' plot_vec(amp_level)]);
                end
                % Give one legend to explain things
                if(p==1)
                    legend(leg_vec);
                end
                ylabel('Amp. Level Prob. ', 'fontsize', 8);   % Do the plot
                H = line([first_on_q(i), first_on_q(i)], [0, 1]); % set the boundary between p and q arms
                set(H, 'LineWidth', 2);  set(H, 'Color', 'r');   % Color it in red
                title(['Pat. : ' HMM_samples{p}], 'fontsize', 8);

                if(p == min(num_samples_in_page, num_samples-(f-1)*num_samples_in_page))  % Temp for 5*5
                    xlabel(['Loc. Chrom. '  num2str(i)] , 'fontsize', 8 );
                end

            end
        end
    end

end


DUMMY=0;


