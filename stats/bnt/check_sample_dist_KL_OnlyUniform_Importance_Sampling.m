% Here we do a very simple thing : Generate a BN.
% Then sample from it many times. Then look at the resulting distribution
% and compute the KL distance between them. This enables us to plot the
% distrubution of the KL for various values of N. 
% Since we want here rare events (KL > eps), we have to use here Importance
% Sampling Technique !!!!!! 

ttttt = cputime;

bnets = 1; % Only one bayesian net so far ...

figure; hold on;
legend_vec = ['ox*+rdgcskm:bv<rdg+cskox*m:bv<v<rdg+cskox*m:bv<v<rdg+cskox*m:bv<'];

max_nsamples = 50000;
res = 2500; % the resulution of N that we use
chunk_size = 2500; % The chunk for saving binomial coefficients 

times_to_sample = 250;
N_samples_vec = res:res:max_nsamples;

nodes_cum_vec = zeros(length(N_samples_vec),40);


% Now if we want to change the number of nodes, and see what we get ...
% We now enumerate on nodes. Note that we use only uniform distribution so
% it is much more efficient, since the joint probability distribution
% decomposes.
for nodes = 1:5;  
    start_nodes = nodes
    
    
    fixed_eps_vec = [0.001, 0.001]; % Now try to evaluate probability to be above some fixed epsilon vector !!!!! 
    number_eps = length(fixed_eps_vec);
        
    epsilon = 0.00000001;
    
    ep = 2.0/3.0;  % probability of an edge in the bnet   
        
    % We now have to adjust Q such that KL(Q||P) = epsilon (approx.)
    % The direction of the Q is random, but for all eps we get the same
    % direction
    
    uni_vec = 0.5 * ones(1, nodes); % representing the uniform distribution
    
    Q_dists = FindIndependentDistWithGivenKLFromUniform(nodes, fixed_eps_vec);  
    Q_cum_dists = cumsum(Q_dists); 
    Q_log_dists = log(Q_dists); Q_one_minus_log_dists = log(1-Q_dists); 
        
%%    samp_dist = zeros( 2^nodes,1); % The sample distribution
    
    % Now sample many times and collect the KL
    KL_Dist_Matrix = zeros(length(N_samples_vec), times_to_sample);   % Gives the sample KL we got in the with N(i) samples in the j-th time 
    Importance_Weight_Matrix_Log = zeros(length(N_samples_vec), times_to_sample);   % Gives the Importance weight P/Q we got in the with N(i) samples in the j-th time
    fixed_eps_upper_cum_matrix = zeros(length(N_samples_vec), number_eps);  % Give the sample probability of the KL to be more than epsilon fro N(i) samples and  fixed_eps(j)
    
    sorted_KL = zeros(6,times_to_sample);
 %%%   all_cum_probs = zeros(2^nodes,6);
    
    % New - straight sampling - no cells !!!! If no memory problems - sample
    % all at once !!!
    %%%alldata = zuk_sample_bnet(bnet, max_nsamples*times_to_sample);
    
%%    samp_dist = zeros( 1,nodes); % The sample distribution
    
    
    % timing checking 
    samptime = 0; proctime = 0;
    
    % Now loop and sample
    for i=1:times_to_sample
        count_dist = zeros( 1,nodes); % The sample distribution
        
        data = rand(nodes,max_nsamples); % sample the i-th node
        
        ind_raw_sample = data > repmat(Q_dists(1,:)', 1, max_nsamples); 
        
        all_raw_sample = zeros(1, max_nsamples); 
        
                
        for cur_node=1:nodes                
            all_raw_sample = all_raw_sample + 2^(cur_node-1) * ind_raw_sample(cur_node,:); %% (data(cur_node,:) > Q_dists(1,cur_node));
        end
        
        
        % Calculate the entropy 
        if(mod(i, 100) == 0)
            do_i = i
        end
        
        % Now we've got the samples. We just need to create the sample distribution
        j=1;
        
        for(cur_samples = N_samples_vec)
            current_raw_sample = all_raw_sample(1:cur_samples);    
            current_sample_hist = hist(current_raw_sample, cur_samples);  % Generate histogram to count how many appear
            current_sample_hist = current_sample_hist(current_sample_hist > 0); % ./ cur_samples; % get rid of all the zeros - Values which never appeared
            
            samptimestop = cputime; samptime = samptime + cputime - samptimestop; samptimestop = cputime; 
            
            % Here calculate the importance weights 
            KL_Dist_Matrix(j,i) = cur_samples * (log(cur_samples)-log(2)*nodes)-sum(current_sample_hist .* log(current_sample_hist)); 
 
            Importance_Weight_Matrix_Log(j, i) = -log(2)*nodes*cur_samples - ...
                sum(sum((1-ind_raw_sample) .* repmat(Q_log_dists(1,:)', 1, max_nsamples))) - sum(sum(ind_raw_sample .* repmat(Q_one_minus_log_dists(1,:)', 1, max_nsamples)));
                        
            % We don't need the constant : - nodes * cur_samples; since we subtruct the min anyway
            
            % Here transfer from counts to distributions
          %%%%  KL_Dist_Matrix(j,i) = -Importance_Weight_Matrix_Log(j, i) ./ cur_samples;  % This is new - use the whole sample
            
            j=j+1;
            
            proctime = proctime+cputime-samptimestop;
        end  % loop on cur_samples 

        KL_Dist_Matrix(:,i) = -KL_Dist_Matrix(:, i) ./ N_samples_vec'; 

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% % %         % Now sample once for the big N to get histogram ..
% % %         uniform_data = rand(nodes,max_nsamples); % sample the i-th node
% % %         uniform_raw_sample = zeros(1, max_nsamples);   
% % %         for cur_node=1:nodes                
% % %             uniform_raw_sample = uniform_raw_sample + 2^(cur_node-1) * (uniform_data(cur_node,:) > 0.5);
% % %         end
% % %         
% % %         uniform_sample_hist = hist(uniform_raw_sample, max_nsamples);  % Generate histogram to count how many appear
% % %         uniform_sample_hist = uniform_sample_hist(uniform_sample_hist > 0); % get rid of all the zeros - Values which never appeared
% % %         
% % %         Impo = 
% % %         
% % %         
% % %         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        
        
    end  % loop on i times to sample 
    
    
    
    % Do it outside to solve the mystery
%     for j=1:length(N_samples_vec)
%         Importance_Weight_Matrix_Log(j, :) = Importance_Weight_Matrix_Log(j, :) -1*log(2)*nodes*N_samples_vec(j);
%     end
    
    % Now update the log weight matrix to avoid underflows
    MinImportanceLogs = min(Importance_Weight_Matrix_Log,[],2); 
    Importance_Weight_Matrix_Log = Importance_Weight_Matrix_Log - repmat(MinImportanceLogs, 1, times_to_sample);
    Importance_Weight_Matrix = exp(Importance_Weight_Matrix_Log); % back to the multiplicative version
    
    % Now try to calculate the prob. to be larger than a GIVEN epsilon 
    for i=1:number_eps
        %%%%%%%%%        fixed_eps_upper_cum_matrix(:,i) = log(mean(KL_Dist_Matrix > fixed_eps_vec(i), 2)) ./ N_samples_vec';  % This is regular (not importance ....)
%         fixed_eps_upper_cum_matrix(:,i) = ...
%             (log(mean(Importance_Weight_Matrix .* (KL_Dist_Matrix > fixed_eps_vec(1)), 2)) + MinImportanceLogs)./ N_samples_vec';  % This is Importance Sampling !!! ....)
        
        % Here do not divide by N, and also add epsilon*N !!!!
        fixed_eps_upper_cum_matrix(:,i) = ...
            (log(mean(Importance_Weight_Matrix .* (KL_Dist_Matrix > fixed_eps_vec(1)), 2)) + MinImportanceLogs) + fixed_eps_vec(1) .* N_samples_vec';  % This is Importance Sampling !!! ....)
        
    end
    

    % Now plot for different N's and for a fixed epsilon
    
    for i=1:1 %%% number_eps
        plot(N_samples_vec, fixed_eps_upper_cum_matrix(:,i), legend_vec(nodes));
    end
    
    nodes_cum_vec(:,nodes) = fixed_eps_upper_cum_matrix(:,i);
    
end % loop on number of nodes 

% Now plot the Sanov Upper Bound
%%%%  plot(N_samples_vec, 2^nodes * log(N_samples_vec+1) ./ N_samples_vec  - fixed_eps_vec(1), '--k');


title('1/N log(Prob. of KL >= eps) for various values of eps and n'); xlabel('N'); ylabel('Prob. KL(N) > eps'); 
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15');


figure; hold on; 
for i=1:length(N_samples_vec)
    plot(nodes_cum_vec(:,i), legend_vec(i)); 
end
title(['prob >= eps as a function of n at eps=' num2str(fixed_eps_vec(1))]); xlabel('n'); ylabel('1/N log prob');

total_simulation_time = cputime - ttttt

total_sampling_time = samptime
total_processing_time = proctime