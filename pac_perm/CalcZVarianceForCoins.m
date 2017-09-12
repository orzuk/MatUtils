% Here we load a data, compute the correlations and get the desired
% fraction
path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
path(path, 'E:\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric\nmm\linalg');


TOL = 0.00000000001;

% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;
rand_flag = GAUSSIAN; 

% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0; 
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0; 
true_corr_flag  = TRUE_AND_SAMPLED;

max_nsamples = 200; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)

res = 20; % the resulution of N that we use
alpha = 0.014; alpha_vec =  [0.012  0.12]; % corresponds to ~70 and ~700 genes

%nsamples=3; 
iters = 10;

nsamples_vec = res:res:max_nsamples;

Ncoins = 1000; Nsamples = size(R.dat, 2) 
res = 1/Ncoins;

coins_p = [res:res:1-res];

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

inf_limit_frac = zeros(2 * length(alpha_vec), length(nsamples_vec));
inf_limit_std = zeros(2 * length(alpha_vec), length(nsamples_vec));
samp_frac_mean = zeros(2 * length(alpha_vec), length(nsamples_vec));
samp_frac_std = zeros(2 * length(alpha_vec), length(nsamples_vec));

index = 1; 
start_plot_ind = plot_ind 


Pearson_corrs = sum(normed_data .* repmat(normed_labels, Ngenes, 1),2);
size(Pearson_corrs)
Fisher_Zs = 0.5 * (log(1+Pearson_corrs) - log(1-Pearson_corrs));


% Now try to fit the best linear distribution to the correlations data.
% We try with the three options to fit : gaussian, uniform and linear.
% So far only Gaussian is supported
if(rand_flag == GAUSSIAN) 
    Fisher_Zs_mean = mean(Fisher_Zs);
    Fisher_Zs_std = std(Fisher_Zs);
    
    % Do fit before the correction !!!! 
    Fisher_Normal_fit_vec = (1/(sqrt(2*pi)*Fisher_Zs_std)) .* exp(-([-0.5:0.001:0.5]-Fisher_Zs_mean).^2 ./ (2.*Fisher_Zs_std^2));
    
    % New ! Make a correction to the std. in order to
    % compensate for the noise of the samples
    Fisher_Zs_std = sqrt(Fisher_Zs_std^2 - (1/(Nsamples-3)));
    
    
    
    
    % Now plot the Fisher correlations and the fit
    figure;  hold on;         
    [hh_vals hh_bins] = hist(Fisher_Zs, 100);
    hist(Fisher_Zs, 100);
    title('Fisher Zs of the data histogram and Gaussian fit');
    hh_bin_size = hh_bins(2)-hh_bins(1);
    %% x_fit_vec = [min(Fisher_Zs):(max(Fisher_Zs)-min(Fisher_Zs))/:max(Fisher_Zs)];
    plot([-0.5:0.001:0.5], Fisher_Normal_fit_vec.*Ngenes.*hh_bin_size, 'r');
    
    % Run over N and calculate the variance of Z
    N_VEC = [10:10:40]; %       N_VEC = [4:floor(Nsamples/2)];
    
    chunk_Pearson_Rhos = zeros(floor(Nsamples/4), Ngenes); % array containing the Zs on chunks
    chunk_Fisher_Zs = zeros(floor(Nsamples/4), Ngenes); % array containing the Zs on chunks
    
    
    
    Z_sig_mean = zeros(1, floor(Nsamples/2));   Z_sig_std = zeros(1, floor(Nsamples/2));
    Rho_sig_mean = zeros(1, floor(Nsamples/2));   Rho_sig_std = zeros(1, floor(Nsamples/2));
    figure; subplot(4,2,1); plot_ind=1;
    for N=N_VEC
        N_IS = N
        chunks_num = floor(Nsamples/N) % number of different blocks
        chunk_Z_sig  = zeros(1,Ngenes);
        chunk_Rho_sig  = zeros(1,Ngenes); 
        
        i=1;
        while(i <= iters)
            curperm = randperm(Nsamples);
            
            update_i=1;
            for j=1:chunks_num
                chunk_normed_labels = R.Labels(curperm((j-1)*N+1:j*N ));
                chunk_normed_labels = chunk_normed_labels-mean(chunk_normed_labels);
                %%%                    - mean(R.Labels(curperm((j-1)*N+1:j*N ))); 
                
                % We need to avoid a case where all the flags 
                % are equal, so we can't define a correlation
                if(min(chunk_normed_labels) == max(chunk_normed_labels))
                    update_i = 0;
                    break;
                end
                
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
                
                if(isempty( find( chunk_Pearson_corrs>=1  )    )  == 0)
                    index = find( chunk_Pearson_corrs>=1  ); index = index(1)
                    chunk_Pearson_corrs(index)
                    data_is = chunk_normed_data(index,:)
                    labels_is  = chunk_normed_labels
                    there_is_one_pearson = 1
                end
                if(isempty(find( chunk_Pearson_corrs<=-1  )    )  == 0)
                    there_is_minus_one_pearson = 1
                end
                %     pear_size = size(chunk_Pearson_corrs)
                %     size(chunk_Fisher_Zs(j,:))
                chunk_Pearson_Rhos(j,:) = chunk_Pearson_corrs;
                chunk_Fisher_Zs(j,:) = 0.5 * (log(1+chunk_Pearson_corrs) - log(1-chunk_Pearson_corrs));
            end
            
            
            %               size(std(chunk_Fisher_Zs(1:chunks_num,:)))
            
            chunk_Z_sig = chunk_Z_sig + std(chunk_Fisher_Zs(1:chunks_num,:));
            chunk_Rho_sig = chunk_Rho_sig + std(chunk_Pearson_Rhos(1:chunks_num,:));
            i_is = i                        
            i=i+update_i;
            
        end
        done_iters = iters
        chunk_Z_sig = chunk_Z_sig./iters;
        Z_sig_mean(N) = mean(chunk_Z_sig);
        Z_sig_std(N) = std(chunk_Z_sig);
        
        chunk_Rho_sig = chunk_Rho_sig./iters;
        Rho_sig_mean(N) = mean(chunk_Rho_sig);
        Rho_sig_std(N) = std(chunk_Rho_sig);
        
          % Plot to see correlation between the 'true' value and the variance
        subplot(4,2,plot_ind); hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. N=' num2str(N)]);
        xlabel('true P'); ylabel('std');
        subplot(4,2,plot_ind+4); plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. N='  num2str(N)]);
        xlabel('true Z'); ylabel('std');

        plot_ind=plot_ind+1;
        %            return
    end
    
    figure;  hold on; title('Z std, from data and from Fisher');
    errorbar(N_VEC, Z_sig_mean(N_VEC),Z_sig_std(N_VEC),'*');
    %%%  plot(N_VEC, Z_sig_mean(N_VEC), '*');                  
    plot(N_VEC, 1./sqrt(N_VEC-3), 'r');  
    xlabel('N'); ylabel('Z std'); legend('from data', 'from Fisher');     
    
    
    figure;  hold on; title('Rho std, from data only');
    errorbar(N_VEC, Rho_sig_mean(N_VEC),Rho_sig_std(N_VEC),'*');
    xlabel('N'); ylabel('Rho std'); %legend('from data', 'from Fisher');     
    
    
    % Plot to see correlation between the 'true' value and the variance
    figure; hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. N=' num2str(N)]);
    xlabel('true P'); ylabel('std');
    figure; hold on; plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. N='  num2str(N)]);
    xlabel('true Z'); ylabel('std');
    
end

% Now calculate the desired N which is needed in order to get some
% desired fraction

% Now plot the results : 


plot_ind = plot_ind+1;

index = index + 1;





