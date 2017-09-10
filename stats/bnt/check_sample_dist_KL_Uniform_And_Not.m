% Here we do a very simple thing : Generate a BN.
% Then sample from it many times. Then look at the resulting distribution
% and compute the KL distance between them. This enables us to plot the
% distrubution of the KL for various values of N. 

ttttt = cputime;

bnets = 1; % Only one bayesian net so far ...
nodes = 4;  % Now if we want to change the number of nodes ... 

max_nsamples = 10000;
res = 200; % the resulution of N that we use

times_to_sample = 250;

fixed_eps_vec = [0.00001, 0.001, 0.01, 0.1]; % Now try to evaluate probability to be above some fixed epsilon vector !!!!! 
number_eps = length(fixed_eps_vec);


epsilon = 0.00000001;

ep = 2.0/3.0;  % probability of an edge in the bnet   

% First generate a random bayesian net
bnet = create_random_bnet(nodes, ep);  % How do we randomize ????? 

% Now try the new sampling method !!!!
bnet_cum_dist = zuk_sample_bnet_prepare(bnet, max_nsamples);


% Now get the distribution from the bnet
% % % engine = enumerative_inf_engine(bnet);
% % % evidence = cell(1, nodes);
% % % engine = enter_evidence(engine, evidence);
% % % 
% % % % Open all the nodes ..
% % % m = marginal_nodes(engine, [1:nodes]);
% % % bnet_dist  = (reshape(m.T,  prod(bnet.node_sizes), 1) + epsilon) / (1.0 + epsilon);

bnet_dist = zuk_bnet_to_probs(bnet);


samp_dist = zeros( 2^nodes,1); % The sample distribution


% Create a random distribution 
ranran = rand(1,16); 
ranran = ranran ./ sum(ranran); 
ranranpermed = ranran(randperm(16));


% Now sample many times and collect the KL
N_samples_vec = res:res:max_nsamples;
KL_Dist_Matrix = zeros(length(N_samples_vec), times_to_sample);   % Gives the sample KL we got in the with N(i) samples in the j-th time 
fixed_eps_upper_cum_matrix = zeros(length(N_samples_vec), number_eps);  % Give the sample probability of the KL to be more than epsilon fro N(i) samples and  fixed_eps(j)



sorted_KL = zeros(6,times_to_sample);
all_cum_probs = zeros(2^nodes,6);

% New - straight sampling - no cells !!!! If no memory problems - sample
% all at once !!!
%%%alldata = zuk_sample_bnet(bnet, max_nsamples*times_to_sample);

for trial = 1:1  % first time bnet, second time, uniform 
samp_dist = zeros( 2^nodes,1); % The sample distribution


% Now loop and sample
for i=1:times_to_sample
   if(trial == 1)
  %%%    data = zuk_sample_bnet(bnet, max_nsamples);  % More time, less memory ...   
      data = zuk_sample_bnet_generate_sample(bnet_cum_dist, max_nsamples);
    
      if(i==1)
          all_cum_probs(:,trial) = bnet_cum_dist; 
      end
  end
  
  if(trial == 2) % Take the same BNET dist ...
        rand_vec = rand(1,max_nsamples);

        data = zeros(1,max_nsamples);

        if(i==1)
            all_cum_probs(:,trial) = bnet_cum_dist;
        end
         % Do by looping on values (not efficient) 
        for val=2^nodes:-1:1
            data(find(rand_vec <= bnet_cum_dist(val))) = val-1; 
        end

  end % if
  
  
  if(trial == 3) % uniform dist ...
         rand_vec = rand(1,max_nsamples);
                 data = zeros(1,max_nsamples);

        bnet_dist = ((1.0/2^nodes) + zeros(1,2^nodes))';
        bnet_cum_dist = [1.0/ 2^nodes : 1.0/ 2^nodes :1]';
        
        if(i==1)
            all_cum_probs(:,trial) = bnet_cum_dist;
        end    
         % Do it by looping on values (not efficient) 
        data = zeros(1,max_nsamples);
        for val=2^nodes:-1:1
            data(find(rand_vec <= bnet_cum_dist(val))) = val-1; 
        end

  end % if

  if(trial == 4) % Fixed non uniform dist
         rand_vec = rand(1,max_nsamples);
        data = zeros(1,max_nsamples);
        bnet_dist = ranran';
        bnet_cum_dist = cumsum(ranran'); %%%%(([1.0/ 2^nodes : 1.0/ 2^nodes :1] .^ 5))';
        
        if(i==1)
            all_cum_probs(:,trial) = bnet_cum_dist;
        end    
         % Do it by looping on values (not efficient) 
        for val=2^nodes:-1:1
            data(find(rand_vec <= bnet_cum_dist(val))) = val-1; 
        end
      
  end
  
   if(trial == 5) % Fixed non uniform dist Again
         rand_vec = rand(1,max_nsamples);
        data = zeros(1,max_nsamples);
        bnet_dist = ranran';
        bnet_cum_dist = cumsum(ranran'); %%%%(([1.0/ 2^nodes : 1.0/ 2^nodes :1] .^ 5))';
        
        if(i==1)
            all_cum_probs(:,trial) = bnet_cum_dist;
        end    
         % Do it by looping on values (not efficient)
        for val=2^nodes:-1:1
            data(find(rand_vec <= bnet_cum_dist(val))) = val-1; 
        end
      
  end

  
  if(trial == 6) % Fixed opposite non uniform dist
  
         rand_vec = rand(1,max_nsamples);
        data = zeros(1,max_nsamples);
         bnet_dist = ranranpermed';
        bnet_cum_dist = (1 -  ([1.0/ 2^nodes : 1.0/ 2^nodes :1] .^ 5))';
        bnet_cum_dist = cumsum(ranranpermed'); %%%%[bnet_cum_dist(end-1:-1:1)', 1]';
        
        if(i==1)
            all_cum_probs(:,trial) = bnet_cum_dist;
        end    
         % Do it by looping on values (not efficient) 
        for val=2^nodes:-1:1
            data(find(rand_vec <= bnet_cum_dist(val))) = val-1; 
        end
      
  end

  
  %%%%%%%  data = alldata((i-1)*max_nsamples+1:i*max_nsamples);
    
    if(mod(i, 100) == 0)
        do_i = i
    end

    
    % Now we've got the samples. We just need to create the sample distribution
    j=1;
    for(cur_samples = (res:res:max_nsamples))
        
        % Work only on part of the data
        cur_data = data(1:cur_samples);   
        
        samp_dist = zeros( 2^nodes,1); % The sample distribution

        % Now we need to transfer the data into a probability distribution
        for val = 0:2^nodes-2
            samp_dist(val+1) = sum(cur_data == val);
        end
        
        % Here we saved a little work for the last one ..
        samp_dist(end) = cur_samples-sum(samp_dist(1:end-1)); 
        
        samp_dist = samp_dist ./ cur_samples;
        
        KL_Dist_Matrix(j,i) = cross_entropy(samp_dist, bnet_dist);
        
        j=j+1;
    end  % loop on cur_samples 
    
end  % loop on i times to sample 


% Now try to calculate the prob. to be larger than a GIVEN epsilon 
if(trial == 1) % Do it only once !!!!
    for i=1:number_eps
        fixed_eps_upper_cum_matrix(:,i) = log(mean(KL_Dist_Matrix > fixed_eps_vec(i), 2)) ./ N_samples_vec';        
    end
end


% Determine the region to look at : 
KL_eps_region = [0:0.000001:0.1];
KL_eps_max = 0.001; % the maximum epsilon we want to display .... 
KL_eps_min = 0.0005;

% Take the rate 1/N times log(Pr(KL > eps))
Sanov_KL_eps_matrix = zeros(length(KL_eps_region), length(N_samples_vec));

colorvec = 'grbmkc';
% Now do the plotting part
for N_ind=1:length(N_samples_vec)
    
%    figure; hold on;  hist( KL_Dist_Matrix(round(end/2),:), 50); title(['KL Histogram for N = '  num2str(N_samples_vec(round(end/2)))]); xlabel('\epsilon'); ylabel('f_{KL}(\epsilon)');

%%%figure; hold on;  hist( KL_Dist_Matrix(N_ind,:), 50); title(['KL Histogram for N = '  num2str(N_samples_vec(N_ind))]); xlabel('\epsilon'); ylabel('f_{KL}(\epsilon)');
        
    % Now plot the cumulative distribution
    sorted_KL(trial,:) = sort(KL_Dist_Matrix(N_ind,:));
    Sanov_KL = log(1-[0:1.0/(length(sorted_KL)-1):1]) ./ max_nsamples;



% % %     % Now we need to compute the rate for the specified region. We 'cut'
% % %     % all that is too big
% % %     cur_KL_end = max(1,sum(sorted_KL < KL_eps_max));
% % %     cur_KL_start = max(1,sum(sorted_KL < KL_eps_min));
% % %     
% % %     sorted_KL = sorted_KL(cur_KL_start:cur_KL_end);
% % %     Sanov_KL =Sanov_KL(cur_KL_start:cur_KL_end);
% % % 
% % % 
% % %     % Now try to determine the scaling law according to Sanov theorem
% % %     plot(sorted_KL, Sanov_KL+0*sorted_KL, colorvec(N_ind)); 
% % %     title(['logprob/N KL 1-Cumulative Distribution for N = '  num2str(N_samples_vec(N_ind))]); xlabel('\epsilon'); ylabel('1/N log (1-F_{KL}(\epsilon))');
% % %     
% % %     
    
    
end % loop on N's

% % % plot(sorted_KL, -sorted_KL, '--');
% % % 
% % % legend([num2str(N_samples_vec) ' lin-fit']);


end % trials loop


  figure; hold on; 
  plot(sorted_KL(1,:), [0:1.0/(length(sorted_KL)-1):1]); 
%  plot(sorted_KL(2,:), [0:1.0/(length(sorted_KL)-1):1], 'm'); 
%  plot(sorted_KL(3,:), [0:1.0/(length(sorted_KL)-1):1], 'r'); 
%  plot(sorted_KL(4,:), [0:1.0/(length(sorted_KL)-1):1], 'g'); 
%  plot(sorted_KL(5,:), [0:1.0/(length(sorted_KL)-1):1], 'c'); 
%  plot(sorted_KL(6,:), [0:1.0/(length(sorted_KL)-1):1], 'k'); 
  title(['KL Cumulative Distribution for N = '  num2str(max_nsamples)]); xlabel('\epsilon'); ylabel('F_{KL}(\epsilon)');
%  legend('Random BNET', 'Same Random BNET', 'Uniform BNET', 'Not Uniform BNET', 'Same Non-Uniform', 'Opposite Non Uniform');
  
  all_cum_probs
  
  
  % Now plot for different N's and for a fixed epsilon
  figure; hold on;
  legend_vec = ['bmrgck'];
  for i=1:number_eps
      plot(N_samples_vec, fixed_eps_upper_cum_matrix(:,i), legend_vec(i));
  end
  title('1/N log(Prob. of KL >= eps) for various values of eps'); xlabel('N'); ylabel('Prob. KL(N) > eps'); 
  legend('0.0001', '0.001', '0.01', '0.1'); 
  
  total_simulation_time = cputime - ttttt