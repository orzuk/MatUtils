% Here we removed most things from the function. We're only interested in
% the average Z variance of each gene, and in the variance of Q as a
% function of N
% The function samples from the data by dividing it to two halfs,
% and computes the empirical distribution (mean/variance/etc.) of the
% Fisher Z score and Pearson P score
function [N_VEC N_HALF_VEC Z_sig_mean Q_sig_mean ] ...
    = CalcTwoZVariancesFromDataFunc(R, alpha_vec, rand_flag, iters, sampling_flag )

global BOOTSTRAP NON_OVERLAP;
global UNIFORM GAUSSIAN LINEAR FROM_DATA mix_GAUSSIAN student_t;
global ONE_SIDE TWO_SIDES;
global TRUE_AND_SAMPLED TWO_SAMPLED;
global TRUE FALSE;

TOL = 0.00000000001;


Ngenes = size(R.dat, 1); Nsamples = size(R.dat, 2)

% Calculate the mean correlation of each gene with survival
normed_labels = R.Labels - mean(R.Labels); normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
normed_data = R.dat - repmat(mean(R.dat, 2), 1, Nsamples);
normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, Nsamples);


corr_mean_vec = sum(repmat(normed_labels, Ngenes, 1) .* normed_data, 2);

% Now try to estimate the variance of the correlation of each gene with
% survival. The method we use is go over all the couples
Fisher_R_to_Z = 1; % Flag saying if to perform fisher RtoZ transformation



Pearson_corrs = sum(normed_data .* repmat(normed_labels, Ngenes, 1),2);
size(Pearson_corrs)
Fisher_Zs = 0.5 * (log(1+Pearson_corrs) - log(1-Pearson_corrs));


% Now try to fit the best linear distribution to the correlations data.
% We try with the three options to fit : gaussian, uniform and linear.
% So far only Gaussian is supported
Fisher_Zs_mean = mean(Fisher_Zs);
Fisher_Zs_std = std(Fisher_Zs);

% Do fit before the correction !!!!
Fisher_Normal_fit_vec = (1/(sqrt(2*pi)*Fisher_Zs_std)) .* exp(-([-0.5:0.001:0.5]-Fisher_Zs_mean).^2 ./ (2.*Fisher_Zs_std^2));

% New ! Make a correction to the std. in order to
% compensate for the noise of the samples
Fisher_Zs_std = sqrt(Fisher_Zs_std^2 - (1/(Nsamples-3)));


% Run over N and calculate the variance of Z
N_VEC = [10:Nsamples]; % Take from five until half just one !!!
N_HALF_VEC = [10:floor(Nsamples/2)];

chunk_Fisher_Zs = zeros(floor(Nsamples/4), Ngenes); % array containing the Zs on chunks



Z_sig_mean = zeros(1, floor(Nsamples/2));   Q_sig_mean = zeros(1,floor(Nsamples/2));


pos_labels = find(R.Labels > 0.5); neg_labels = find(R.Labels < 0.5);
pos_num = length(pos_labels); neg_num = Nsamples-pos_num; % Negative and positive labels
[val_labels ind_labels] = sort(R.Labels);

for N=N_VEC
    N_IS = N
    chunks_num = floor(Nsamples/N); % number of different blocks
    chunk_Z_sig  = zeros(1,Ngenes);

    Q_chunk_sig = 0; 
    
    i=1;

    % New ! Do not allow more than two chunks !!!! 
    chunks_num=min(chunks_num,2); 
    
    % randomize a subset many times and each time measure the overlap
    while(i <= iters)

        if(mod(i,100) == 1)
            cur_iter = i
        end

        
        if(sampling_flag == BOOTSTRAP) % Here sample with replacement
            curperm = ceil(Nsamples*rand(1,Nsamples));  
            
            % Now make sure that we have at least one of each set            
            curperm(1:N:1+N*(chunks_num-1)) = pos_labels(ceil(pos_num*rand(1,chunks_num)));
            curperm(2:N:2+N*(chunks_num-1)) = neg_labels(ceil(neg_num*rand(1,chunks_num)));
            
        else
            posperm = randperm(pos_num);  negperm = randperm(neg_num);
            curperm = ind_labels([negperm posperm + neg_num]);



            order_perm = zeros(1,Nsamples);
            order_perm(1:N:1+N*(chunks_num-1)) = 1:chunks_num;
            order_perm(2:N:2+N*(chunks_num-1)) = Nsamples:-1:Nsamples-chunks_num + 1;

            small_order_perm = randperm(Nsamples-2*chunks_num) + chunks_num;
            bad_indexes = [1:N:1+N*(chunks_num-1) 2:N:2+N*(chunks_num-1)];
            small_ind_set = setdiff(1:Nsamples, bad_indexes);
            order_perm(small_ind_set) = small_order_perm;

            curperm = curperm(order_perm);
        end

        % Reorganize the labels. Now we know that the ones are at the end
       % temp_labels = R.Labels(curperm)  
        
        
            
        update_i=1;

        % Do preliminary check before everything - is still needed?????
% % % % % % % % %         for j=1:chunks_num
% % % % % % % % %             chunk_normed_labels = R.Labels(curperm((j-1)*N+1:j*N ));
% % % % % % % % %             chunk_normed_data = R.dat(:,curperm((j-1)*N+1:j*N ));
% % % % % % % % % 
% % % % % % % % %             % We need to avoid a case where all the flags or the
% % % % % % % % %             % data's!!! are equal, so we can't define a correlation
% % % % % % % % %             if(min(chunk_normed_labels) == max(chunk_normed_labels))
% % % % % % % % %                 update_i = 0;
% % % % % % % % %                 bad_chunk_is = j
% % % % % % % % %                 chunks_num_is = chunks_num
% % % % % % % % %                 problem_is_in_labels = 1
% % % % % % % % %                 break;
% % % % % % % % %             end
% % % % % % % % %             if(~isempty(find(min(chunk_normed_data,[],2) == max(chunk_normed_data,[],2)))) % We simultaniously check for all rows
% % % % % % % % %                 update_i = 0;
% % % % % % % % %                 problem_is_in_data = 199
% % % % % % % % %                 indexes_f = find(min(chunk_normed_data,[],2) == max(chunk_normed_data,[],2));
% % % % % % % % %                 min_f = min(chunk_normed_data,[],2);
% % % % % % % % %                 BAD_VALS = min_f(indexes_f)
% % % % % % % % %                 break;
% % % % % % % % %             end
% % % % % % % % % 
% % % % % % % % %         end

        % Getting here means that the vectors are not constant !
        if(update_i == 1)
            % Now do real loop
            for j=1:chunks_num
                chunk_normed_labels = R.Labels(curperm((j-1)*N+1:j*N ));
                chunk_normed_labels = chunk_normed_labels-mean(chunk_normed_labels);
                chunk_normed_labels = chunk_normed_labels ./ sqrt(sum(chunk_normed_labels.^2));
                chunk_normed_data = R.dat(:,curperm((j-1)*N+1:j*N ));
                chunk_normed_data = chunk_normed_data - repmat(mean(chunk_normed_data, 2), 1, N);
                chunk_normed_data = chunk_normed_data ./ repmat(sqrt(sum(chunk_normed_data.^2, 2)),1, N);
                chunk_Pearson_corrs = (sum(chunk_normed_data .* repmat(chunk_normed_labels, Ngenes, 1),2))';
                chunk_Pearson_corrs = min(chunk_Pearson_corrs,1-TOL);  chunk_Pearson_corrs = max(chunk_Pearson_corrs,-1+TOL);

                
%                chunk_Fisher_Zs(j,:) = 0.5 * (log(1+chunk_Pearson_corrs) - log(1-chunk_Pearson_corrs));  % Fisher-Z-transform
                chunk_Fisher_Zs(j,:) = atanh(chunk_Pearson_corrs);
            end

            % Get the variance from the measurements
            if(chunks_num > 1)  % This is only up to Nsamples/2
                chunk_Z_sig = chunk_Z_sig + var(chunk_Fisher_Zs(1:chunks_num,:),0,1);
            end

            Q_chunk_sig = Q_chunk_sig + mean(var(chunk_Fisher_Zs(1:chunks_num,:),0,2));
            

            if(mod(i,50) == 0)
                i_is = i
            end
            i=i+update_i;

        end  % if update i

    end
    done_iters = iters
    
    if(chunks_num > 1)
        chunk_Z_sig = sqrt(chunk_Z_sig./iters); % Now move from variance to std.
        Z_sig_mean(N) = sqrt(mean(chunk_Z_sig)); % Do sqrt after mean!!!! 
    end
    
    Q_sig_mean(N) = sqrt(Q_chunk_sig/iters); % Now move from variance to std.
    
end

% Take only the relevant values
Z_sig_mean = Z_sig_mean(N_HALF_VEC);
Q_sig_mean = Q_sig_mean(N_VEC);






