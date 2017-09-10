% The function samples from the data by dividing it to two halfs,
% and computes the empirical distribution (mean/variance/etc.) of the
% Fisher Z score and Pearson P score
function [N_VEC Pearson_corrs Fisher_Zs Rho_sig_mean Rho_sig_std Rho_bias_mean Rho_bias_std ...
    Z_sig_mean Z_sig_std Z_bias_mean Z_bias_std chunk_Rho_sig chunk_Z_sig chunk_Rho_bias chunk_Z_bias ...
    hist_Pearson_Ps hist_Fisher_Zs random_genes_picked] ...
    = CalcZVarianceFromDataFunc(R, alpha_vec, rand_flag, iters )

TOL = 0.00000000001;

UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; FROM_DATA = 3; mix_GAUSSIAN = 4;
% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0; 


% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0; 


Ngenes = size(R.dat, 1); Nsamples = size(R.dat, 2) 

% Calculate the mean correlation of each gene with survival
normed_labels = R.Labels - mean(R.Labels); normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
normed_data = R.dat - repmat(mean(R.dat, 2), 1, Nsamples); 
normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, Nsamples);
%normed_data = normed_data';


corr_mean_vec = sum(repmat(normed_labels, Ngenes, 1) .* normed_data, 2); 

% Now try to estimate the variance of the correlation of each gene with
% survival. The method we use is go over all the couples 
Fisher_R_to_Z = 1; % Flag saying if to perform fisher RtoZ transformation

plot_ind=1;

index = 1; 
start_plot_ind = plot_ind 


Pearson_corrs = sum(normed_data .* repmat(normed_labels, Ngenes, 1),2);
size(Pearson_corrs)
Fisher_Zs = 0.5 * (log(1+Pearson_corrs) - log(1-Pearson_corrs));


% Now try to fit the best linear distribution to the correlations data.
% We try with the three options to fit : gaussian, uniform and linear.
% So far only Gaussian is supported
%%if(rand_flag == GAUSSIAN) 
    Fisher_Zs_mean = mean(Fisher_Zs);
    Fisher_Zs_std = std(Fisher_Zs);
    
    % Do fit before the correction !!!! 
    Fisher_Normal_fit_vec = (1/(sqrt(2*pi)*Fisher_Zs_std)) .* exp(-([-0.5:0.001:0.5]-Fisher_Zs_mean).^2 ./ (2.*Fisher_Zs_std^2));
    
    % New ! Make a correction to the std. in order to
    % compensate for the noise of the samples
    Fisher_Zs_std = sqrt(Fisher_Zs_std^2 - (1/(Nsamples-3)));
    
    
    
    
    % Now plot the Fisher correlations and the fit
% % % % % %     figure;  hold on;         
% % % % % %     [hh_vals hh_bins] = hist(Fisher_Zs, 100);
% % % % % %     hist(Fisher_Zs, 100);
% % % % % %     title('Fisher Zs of the data histogram and Gaussian fit');
% % % % % %     hh_bin_size = hh_bins(2)-hh_bins(1);
% % % % % %     %% x_fit_vec = [min(Fisher_Zs):(max(Fisher_Zs)-min(Fisher_Zs))/:max(Fisher_Zs)];
% % % % % %     plot([-0.5:0.001:0.5], Fisher_Normal_fit_vec.*Ngenes.*hh_bin_size, 'r');
    
    % Run over N and calculate the variance of Z
    N_VEC = [floor(Nsamples/2):floor(Nsamples/2)] % Take just one !!! 
    
    chunk_Pearson_Rhos = zeros(floor(Nsamples/4), Ngenes); % array containing the Zs on chunks
    chunk_Fisher_Zs = zeros(floor(Nsamples/4), Ngenes); % array containing the Zs on chunks
    
    
    
    Z_sig_mean = zeros(1, floor(Nsamples/2));   Z_sig_std = zeros(1, floor(Nsamples/2));
    Rho_sig_mean = zeros(1, floor(Nsamples/2));   Rho_sig_std = zeros(1, floor(Nsamples/2));
    Z_bias_mean = zeros(1, floor(Nsamples/2));   Z_bias_std = zeros(1, floor(Nsamples/2));
    Rho_bias_mean = zeros(1, floor(Nsamples/2));   Rho_bias_std = zeros(1, floor(Nsamples/2));
    
    %% figure; subplot(4,2,1); 
    plot_ind=1;
    for N=N_VEC
        N_IS = N
        chunks_num = floor(Nsamples/N) % number of different blocks
        chunk_Z_sig  = zeros(1,Ngenes);
        chunk_Rho_sig  = zeros(1,Ngenes);
        chunk_Z_bias  = zeros(1,Ngenes);
        chunk_Rho_bias  = zeros(1,Ngenes);

        i=1;
        
        % We also collect ALL the obtained correlations!!! 
        NUM_RAND_GENES = 36;
        hist_Pearson_Ps = zeros(NUM_RAND_GENES,iters);
        hist_Fisher_Zs = zeros(NUM_RAND_GENES,iters);
        
        % Pick 16 genes at random forwhich we calculate the entire
        % distribution
        random_genes_picked = randperm(Ngenes); random_genes_picked = random_genes_picked(1:NUM_RAND_GENES);
        
        
        % randomize a subset many times and each time measure the overlap
        while(i <= iters)

            if(mod(i,100) == 1)
                cur_iter = i
            end
            curperm = randperm(Nsamples);

            update_i=1;

            % Do preliminary check before everything
            for j=1:chunks_num
                chunk_normed_labels = R.Labels(curperm((j-1)*N+1:j*N ));
                chunk_normed_data = R.dat(:,curperm((j-1)*N+1:j*N ));

                
                % We need to avoid a case where all the flags or the
                % data's!!! are equal, so we can't define a correlation
                if(min(chunk_normed_labels) == max(chunk_normed_labels))
                    update_i = 0;
                    problem_is_in_labels = 1
                    break;
                end
                if(~isempty(find(min(chunk_normed_data,[],2) == max(chunk_normed_data,[],2)))) % We simultaniously check for all rows                  
                    update_i = 0;
                    %%% We don't need to print now. Just move on
                  %%%  problem_is_in_data = 199
                    indexes_f = find(min(chunk_normed_data,[],2) == max(chunk_normed_data,[],2));
                    min_f = min(chunk_normed_data,[],2);
                  %%%  BAD_VALS = min_f(indexes_f)
                    break;
                end
                
            end

            % Getting here means that the vectors are not constant !
            if(update_i == 1)
                % Now do real loop
                for j=1:chunks_num
                    chunk_normed_labels = R.Labels(curperm((j-1)*N+1:j*N ));
                    chunk_normed_labels = chunk_normed_labels-mean(chunk_normed_labels);
                    %%%                    - mean(R.Labels(curperm((j-1)*N+1:j*N )));


                    %                     chunk_normed_labels
                    %                     sqrt(sum(chunk_normed_labels.^2))
                    chunk_normed_labels = chunk_normed_labels ./ sqrt(sum(chunk_normed_labels.^2));
                    %    size_R_dat = size(R.dat)
                    %    size(R.dat(:,curperm((j-1)*N+1:j*N)))
                    %    size(repmat(mean(R.dat(:,curperm((j-1)*N+1:j*N )), 2), 1, N))
                    chunk_normed_data = R.dat(:,curperm((j-1)*N+1:j*N ));
                    chunk_normed_data = chunk_normed_data - repmat(mean(chunk_normed_data, 2), 1, N);
                    chunk_normed_data = chunk_normed_data ./ repmat(sqrt(sum(chunk_normed_data.^2, 2)),1, N);
                    chunk_Pearson_corrs = (sum(chunk_normed_data .* repmat(chunk_normed_labels, Ngenes, 1),2))';
                    chunk_Pearson_corrs = min(chunk_Pearson_corrs,1-TOL);  chunk_Pearson_corrs = max(chunk_Pearson_corrs,-1+TOL);

               
                    
                    % Remove for now for slightly better performance 
% % % % % % % % % % % % %                     if(isempty( find( chunk_Pearson_corrs>=1  )    )  == 0)
% % % % % % % % % % % % %                         index = find( chunk_Pearson_corrs>=1  ); index = index(1)
% % % % % % % % % % % % %                         chunk_Pearson_corrs(index)
% % % % % % % % % % % % %                         data_is = chunk_normed_data(index,:)
% % % % % % % % % % % % %                         labels_is  = chunk_normed_labels
% % % % % % % % % % % % %                         there_is_one_pearson = 1
% % % % % % % % % % % % %                     end
% % % % % % % % % % % % %                     if(isempty(find( chunk_Pearson_corrs<=-1  )    )  == 0)
% % % % % % % % % % % % %                         there_is_minus_one_pearson = 1
% % % % % % % % % % % % %                     end



%     pear_size = size(chunk_Pearson_corrs)
                    %     size(chunk_Fisher_Zs(j,:))
                    chunk_Pearson_Rhos(j,:) = chunk_Pearson_corrs;
                    chunk_Fisher_Zs(j,:) = 0.5 * (log(1+chunk_Pearson_corrs) - log(1-chunk_Pearson_corrs));  % Fisher-Z-transform
                end

                % Use only the last chunk - quite wastefull
                hist_Pearson_Ps(:,i) = chunk_Pearson_corrs(random_genes_picked);
                hist_Fisher_Zs(:,i) = chunk_Fisher_Zs(j,random_genes_picked);
                
                
                %               size(std(chunk_Fisher_Zs(1:chunks_num,:)))

                % Get the variance from the measurements
                chunk_Z_sig = chunk_Z_sig + std(chunk_Fisher_Zs(1:chunks_num,:),0, 1) .^ 2;
                chunk_Rho_sig = chunk_Rho_sig + std(chunk_Pearson_Rhos(1:chunks_num,:),0, 1) .^ 2;


                % Get the bias from the measurements
                chunk_Z_bias = chunk_Z_bias + mean(chunk_Fisher_Zs(1:chunks_num,:));
                chunk_Rho_bias = chunk_Rho_bias + mean(chunk_Pearson_Rhos(1:chunks_num,:));

                if(mod(i,50) == 0)
                    i_is = i
                end
                i=i+update_i;

            end  % if update i

        end
        done_iters = iters
        chunk_Z_sig = sqrt(chunk_Z_sig./iters);
        Z_sig_mean(N) = mean(chunk_Z_sig);
        Z_sig_std(N) = std(chunk_Z_sig);

        chunk_Rho_sig = sqrt(chunk_Rho_sig./iters);
        Rho_sig_mean(N) = mean(chunk_Rho_sig);
        Rho_sig_std(N) = std(chunk_Rho_sig);

        chunk_Rho_bias = chunk_Rho_bias./iters - Pearson_corrs';
        Rho_bias_mean(N) = mean(chunk_Rho_bias);
        Rho_bias_std(N) = std(chunk_Rho_bias);

        chunk_Z_bias = chunk_Z_bias./iters - Fisher_Zs';
        Z_bias_mean(N) = mean(chunk_Z_bias);
        Z_bias_std(N) = std(chunk_Z_bias);


        % Plot to see correlation between the 'true' value and the variance
        %%        subplot(4,2,plot_ind); hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. N=' num2str(N)]);
        %%        xlabel('true P'); ylabel('std');
        %%       subplot(4,2,plot_ind+4); plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. N='  num2str(N)]);
        %%      xlabel('true Z'); ylabel('std');

        plot_ind=plot_ind+1;
        %            return
    end

    
%%end % if Gaussian

% Now calculate the desired N which is needed in order to get some
% desired fraction

% Now plot the results : 


plot_ind = plot_ind+1;

index = index + 1;


% % Now plot the histograms !!!
% figure; subplot(4,4,1);
% 
% for i=1:4
%     for j=1:4
%         subplot(4,4,4*(i-1)+j); hold on;
%         hist(hist_Pearson_Ps(4*(i-1)+j,:), 100); title(['Pearson hist gene # ' num2str(random_genes_picked(4*(i-1)+j))]);
%     end
% end
% 
% 
% figure; subplot(4,4,1);
% 
% for i=1:4
%     for j=1:4
%         subplot(4,4,4*(i-1)+j); hold on;
%         hist(hist_Fisher_Zs(4*(i-1)+j,:), 100); title(['Fisher hist gene # ' num2str(random_genes_picked(4*(i-1)+j))]);
%     end
% end





