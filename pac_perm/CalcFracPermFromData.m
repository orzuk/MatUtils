% Here we load a data, compute the correlations and get the desired
% fraction
path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
path(path, 'E:\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric\nmm\linalg');


KNN = 10; % parameter for knn missing values algorithm
Fisher_R_to_Z = 1; % Flag saying if to perform fisher RtoZ transformation

% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;
rand_flag = GAUSSIAN;

% data files
OLD_VANT_VEER = 1; NEW_ROSSETA = 2; WANG = 3;     % gene expression
TOPIC = 4; % Author-topic matching
NIPS_DOROTHEA = 5; NIPS_ARCENE = 6; NIPS_DEXTER = 7; NIPS_GISETTE = 8; % Various datasetes from NIPS contest

% Vector saying which datasets come in a sparse form
IS_SPARSE_VEC = [0,0,0,0,0,0,0,0];

TRUE = 1; FALSE = 0;
calc_corrs_flag = FALSE; % flag saying if to calculate the correlations from the data, or do something else


% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TRUE_AND_SAMPLED;

max_nsamples = 100; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)

res = 25; % the resulution of N that we use
alpha = 0.014; alpha_vec =  [0.012  0.12]; % corresponds to ~70 and ~700 genes
TOP_GENES = 70;


%nsamples=3;
iters = 100;

nsamples_vec = res:res:max_nsamples;

% Now chose which data to load - New : Run on all data's !!!!


for true_corr_flag=0:1

    plot_ind=1;
    figure; subplot(2,2,1); hold on;

    for data_flag = TOPIC % NIPS_GISETTE; %NIPS_DEXTER; %%NIPS_ARCENE; %NIPS_DOROTHEA; % TOPIC; % WANG; %%% NEW_ROSSETA; %%%OLD_VANT_VEER;

        clear R; % free memory


        if(data_flag == OLD_VANT_VEER)
            load('breastData.mat'); % The name of the data file containing the expression and labels.
            data_str = 'VantVeerOld';
        end
        if(data_flag == NEW_ROSSETA)
            loading_big_rosseta = 1
            load('Rosetta_data.mat'); % The name of the data file containing the expression and labels.
            data_str = 'VantVeerNew';

            finished_loading_big_rosseta = 1
            R.dat = knnimpute(data, KNN);
            R.Labels = cell2mat(samples); R.Labels = R.Labels';
            loading_big_rosseta = 1
        end

        if(data_flag == WANG)
            load('WANG_DATA.mat');
            data_str = 'Wang';

            % Here we must do log since Gary didn't do it
            R.dat = log(gene_expression);
            R.Labels = Prognosis';
        end

        if(data_flag == TOPIC)
            load('topic_data.mat');
            data_str = 'Author-Topic';
            R.dat = Topics_N_counts;
            R.Labels = N_Labels';
            Pearson_corrs = words_corr;   % Many zero correlations in this dataset.
        end

        if(data_flag == NIPS_DOROTHEA)
            load('an_dorothea.mat');
            data_str = 'NIPS-Dorothea';
            R.dat = T;
            R.Labels = Labels';
            Pearson_corrs = corr;   % Note : Here these aren't correlations but something else ! Distribution not gaussian but a lot of mass on the right side
        end

        if(data_flag == NIPS_ARCENE)
            load('ARCENE_data_file.mat');
            data_str = 'NIPS-Arcene';
            R.dat = T;
            R.Labels = Labels';
            Pearson_corrs = corr;   % Very weak correlations in this dataset
        end

        if(data_flag == NIPS_DEXTER)
            load('DEXTER_data_file.mat');
            data_str = 'NIPS-Dexter';
            R.dat = T;
            R.Labels = Labels';
            Pearson_corrs = corr;   % Very weak correlations in this dataset. Most are zero. NOT a gaussian distribution
        end

        if(data_flag == NIPS_GISETTE)
            load('GISETTE_data_file.mat');
            data_str = 'NIPS-Gisette';
            R.dat = T;
            R.Labels = Labels';
            Pearson_corrs = corr;   % Very weak correlations in this dataset. Not symmetric, mostly negative correlations.
        end



        Ngenes = size(R.dat, 1)
        Nsamples = size(R.dat, 2)


        if(calc_corrs_flag == TRUE)
            % Calculate the mean correlation of each gene with survival
            normed_labels = R.Labels - mean(R.Labels); normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
            normed_data = R.dat - repmat(mean(R.dat, 2), 1, Nsamples);
            normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, Nsamples);
            %normed_data = normed_data';

            Pearson_corrs = sum(normed_data .* repmat(normed_labels, Ngenes, 1),2);
            size(Pearson_corrs)

        end % here the correlations are already given


        % Now try to estimate the variance of the correlation of each gene with
        % survival. The method we use is go over all the couples

        inf_limit_frac = zeros(2 * length(alpha_vec), length(nsamples_vec));
        inf_limit_std = zeros(2 * length(alpha_vec), length(nsamples_vec));
        samp_frac_mean = zeros(2 * length(alpha_vec), length(nsamples_vec));
        samp_frac_std = zeros(2 * length(alpha_vec), length(nsamples_vec));

        index = 1;

        %%%     for beta = alpha_vec % alpha_vec alpha
        alpha = alpha_vec(1); % TOP_GENES / Ngenes; % make the number of genes always constant (70)
        true_flag_is = true_corr_flag
        alpha_is = alpha
        start_plot_ind = plot_ind

        if(Fisher_R_to_Z)


            Fisher_Zs = 0.5 * (log(1+Pearson_corrs) - log(1-Pearson_corrs));


            SAVE_FISHER_Z{data_flag+1} = Fisher_Zs; % Save the Fisher Z's for plotting later

            % Now try to fit the best linear/gaussian/whatever distribution to the correlations data.
            % We try with the three options to fit : gaussian, uniform and linear.
            % So far only Gaussian is supported
            if(rand_flag == GAUSSIAN)
                Fisher_Zs_mean = mean(Fisher_Zs);
                Fisher_Zs_std = std(Fisher_Zs);

                % New ! Make a correction to the std. in order to
                % compensate for the noise of the samples
                Fisher_Zs_std = sqrt(Fisher_Zs_std^2 - (1/(Nsamples-3)));

                Fisher_Normal_fit_vec = (1/(sqrt(2*pi)*Fisher_Zs_std)) .* exp(-([-0.5:0.001:0.5]-Fisher_Zs_mean).^2 ./ (Fisher_Zs_std^2));


            end

            % Now calculate the desired N which is needed in order to get some
            % desired fraction
            i=1;
            for(nsamples=nsamples_vec)

                done_nsamples = nsamples

                % Compute analytically the mean and (approx. ) variance
                [inf_limit_frac(index, i) inf_limit_std(index, i) x_alpha ] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, nsamples-3, alpha);  % We use here the approximation of Fisher



                % Now do also sampling and compare to analytic result
                kept_frac_dist = sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, nsamples-3, Ngenes, alpha, iters, R, FALSE);
                samp_frac_mean(index, i) = mean(kept_frac_dist)/(alpha*Ngenes); samp_frac_std(index, i) = std(kept_frac_dist./(alpha*Ngenes));
                i=i+1;

            end

            % adjust for Ngenes
            inf_limit_std(index,:) = inf_limit_std(index,:) ./ sqrt(Ngenes);


            % Now plot the results :
            % Now print the figure. The simulation together with the delta of ngenes->infinity
            if(one_side_flag)
                side_str = 'one sided, ';
            else
                side_str = 'two sides, ';
            end
            if(true_corr_flag==TRUE_AND_SAMPLED)
                corr_str = 'true&samp., ';
            else
                corr_str = 'two samp., ';
            end
            if(rand_flag == GAUSSIAN)
                rand_str = 'gauss., ';
            end
            if(rand_flag == UNIFORM)
                rand_str = 'uni., ';
            end
            if(rand_flag == LINEAR)
                rand_str = 'lin., ';
            end




            subplot(2,2,plot_ind); hold on;
            if(true_corr_flag == TRUE_AND_SAMPLED)
                errorbar(nsamples_vec, inf_limit_frac(index,:), inf_limit_std(index,:));
            else
                plot(nsamples_vec, inf_limit_frac(index,:));
            end
            errorbar(nsamples_vec, samp_frac_mean(index,:), samp_frac_std(index,:), 'r');   ylabel('Mean kept Frac f*'); xlabel('No. Samples'); legend( 'Saddle approx.','sampled', 2);
            title([data_str ' Data, f dist. ' rand_str corr_str side_str 'Ngenes=' num2str(Ngenes) ' alpha=' num2str(alpha) ' iters=' num2str(iters)]);


        else  % Here do some sampling-based method


            corr_num_iters = 10000; % number of iterations to randomize couples


            for iter = 1:corr_num_iters
                % randomize the couples
                rand_perm = randperm(1,Nsamples);

                for i=1:2:floor(Nsamples/2)*2-1
                    % Calculate the current two-point data and normalize it
                    cur_data_vec = R.data([rand_perm(i), rand_perm(i+1)],:);
                    cur_data_vec = cur_data_vec - repmat(mean(cur_data_vec, 2), 1, 2);
                    cur_data_vec = cur_data_vec ./ repmat( sqrt(sum( cur_data_vec .^2, 1)), 1, 2);

                    cur_labels_vec = R.labels(rand_perm(i), rand_perm(i+1));

                    % Calculate the correlations
                    curr_corr_vec = cur_data_vec .* repmat(cur_labels_vec, 1, Ngenes);

                    for j=1:Ngenes % for now do correlation gene by gene !

                        garbage  = curr_corr_vec;

                    end

                end
            end
        end % if Fisher flag

        plot_ind = plot_ind+1;

        index = index + 1;

        %%%      end  % loop on alpha
        %%% end  % loop on two samps flag


    end  % end loop on data types


end % loop on true corrs



% Now plot the Fisher Z's for comparison
figure; subplot(2,2,1); hold on;
for i=1:3
    data_flag = i-1;

    if(data_flag == OLD_VANT_VEER)
        data_str = 'VantVeerOld';
    end
    if(data_flag == NEW_ROSSETA)
        data_str = 'VantVeerNew';
    end
    if(data_flag == WANG)
        data_str = 'Wang';
    end
    subplot(2,2,i); hist(SAVE_FISHER_Z{i}, 250); xlabel('Correlation'); ylabel('Freq.'); title([data_str ' Data Hist. of all genes Fisher Zs']);
end














return; % For now dont do std.s .


% Now plot only the standard deviations
figure; subplot(2,2,1); hold on;
index = 1;
for true_corr_flag=1:1 % 0:1
    for alpha = alpha_vec % alpha_vec


        % Now print the figure. The simulation together with the delta of ngenes->infinity
        if(one_side_flag)
            side_str = 'one sided, ';
        else
            side_str = 'two sides, ';
        end
        if(true_corr_flag==TRUE_AND_SAMPLED)
            corr_str = 'true&sampled, ';
        else
            corr_str = 'two sampled, ';
        end
        if(rand_flag == GAUSSIAN)
            rand_str = 'gaussian, ';
        end
        if(rand_flag == UNIFORM)
            rand_str = 'uniform, ';
        end
        if(rand_flag == LINEAR)
            rand_str = 'linear, ';
        end

        if(data_flag == OLD_VANT_VEER)
            data_str = 'VantVeerOld';
        end
        if(data_flag == NEW_ROSSETA)
            data_str = 'VantVeerNew';
        end
        if(data_flag == WANG)
            data_str = 'Wang';
        end

        subplot(2,2,index); hold on;

        errorbar(nsamples_vec, inf_limit_frac(2+index,:), inf_limit_std(2+index,:));
        errorbar(nsamples_vec, samp_frac_mean(2+index,:), samp_frac_std(2+index,:), 'r');   ylabel('Mean kept Frac f*'); xlabel('No. Samples'); legend( 'Saddle approx.','sampled', 2);
        title([data_str ' Data, f dist. ' rand_str corr_str side_str 'Ngenes=' num2str(Ngenes) ' alpha=' num2str(alpha) ' iters=' num2str(iters)]);


        %        plot(nsamples_vec, inf_limit_std(index,:));
        %        plot(nsamples_vec, samp_frac_std(index,:), 'r');   ylabel('Std kept Frac f*'); xlabel('No. Samples'); legend( 'saddle approx.','sampled', 2);
        %        title(['Kept frac. Std. ' rand_str corr_str side_str 'Ngenes=' num2str(Ngenes) ' alpha=' num2str(alpha) ' iters=' num2str(iters)]);

        index = index +1;

    end
end






return;


% Now show as an example the bahaviour for one nsamples value
iters = 5000; true_corr_flag = TRUE_AND_SAMPLED; nsamples = 300; alpha  = alpha_vec(1); f_res = 0.0005; soft_constrain_flag = 0;
alpha
CheckPermTopFunc(rand_flag, one_side_flag, true_corr_flag, Fisher_Zs_std, nsamples-3, Ngenes, alpha, iters, f_res, soft_constrain_flag, miu, prior);

% Where is the plot ???? 

