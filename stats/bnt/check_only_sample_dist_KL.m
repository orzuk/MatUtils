% Here we do a very simple thing : Generate a BN.
% Then sample from it many times. Then look at the resulting distribution
% and compute the KL distance between them. This enables us to plot the
% distrubution of the KL for various values of N. 
bnets = 1; % Only one bayesian net so far ...
nodes = 4; 

max_nsamples = 1000;
res = 500;

times_to_sample = 5000;


epsilon = 0.00000001;

ep = 2.0/3.0;  % probability of an edge in the bnet   

% First generate a random bayesian net
bnet = create_random_bnet(nodes, ep);  % How do we randomize ????? 

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



% Now sample many times and collect the KL
N_samples_vec = res:res:max_nsamples;
KL_Dist_Matrix = zeros(length(N_samples_vec), times_to_sample);





% New - straight sampling - no cells !!!! If no memory problems - sample
% all at once !!!
%%%alldata = zuk_sample_bnet(bnet, max_nsamples*times_to_sample);

% Now loop and sample
for i=1:times_to_sample
 
  data = zuk_sample_bnet(bnet, max_nsamples);  % More time, less memory ...   
  %%%%%%%  data = alldata((i-1)*max_nsamples+1:i*max_nsamples);
    
    if(mod(i, 100) == 0)
        do_i = i
    end

% % %     % Old slow sampling ... 
% % %     sam = cell(nodes, max_nsamples);
% % %     
% % %     % sample data from it. For now do it one by one ! (very slow ...) 
% % %     for j = 1:max_nsamples
% % %         sam(:,j) = sample_bnet(bnet);
% % %     end
% % %     
% % %     % transfer from cell to a number matrix
% % %     data = cell2num(sam);
% % %     
% % %     % Transfer the data to numbers between 0 to 2^nodes-1
% % %     data = data-1;
% % %     power_vec = (repmat(2.^[0:nodes-1], max_nsamples, 1))';
% % %     data = sum(data .* power_vec);
    

    
    % Now we've got the samples. We just need to create the sample distribution
    j=1;
    for(cur_samples = (res:res:max_nsamples))
        
        % Work only on part of the data
        cur_data = data(1:cur_samples);   
        
        
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
    
end


% Determine the region to look at : 
KL_eps_region = [0:0.000001:0.1];
KL_eps_max = 0.001; % the maximum epsilon we want to display .... 
KL_eps_min = 0.0005;

% Take the rate 1/N times log(Pr(KL > eps))
Sanov_KL_eps_matrix = zeros(length(KL_eps_region), length(N_samples_vec));

colorvec = 'grbmkc';
figure; hold on;
% Now do the plotting part
for N_ind=1:length(N_samples_vec)
    
%    figure; hold on;  hist( KL_Dist_Matrix(round(end/2),:), 50); title(['KL Histogram for N = '  num2str(N_samples_vec(round(end/2)))]); xlabel('\epsilon'); ylabel('f_{KL}(\epsilon)');

%%%figure; hold on;  hist( KL_Dist_Matrix(N_ind,:), 50); title(['KL Histogram for N = '  num2str(N_samples_vec(N_ind))]); xlabel('\epsilon'); ylabel('f_{KL}(\epsilon)');
        
    % Now plot the cumulative distribution
    sorted_KL = sort(KL_Dist_Matrix(N_ind,:));
    Sanov_KL = log(1-[0:1.0/(length(sorted_KL)-1):1]) ./ max_nsamples;

  figure; hold on; plot(sorted_KL, [0:1.0/(length(sorted_KL)-1):1]); title(['KL Cumulative Distribution for N = '  num2str(max_nsamples)]); xlabel('\epsilon'); ylabel('F_{KL}(\epsilon)');

    % Now we need to compute the rate for the specified region. We 'cut'
    % all that is too big
    cur_KL_end = max(1,sum(sorted_KL < KL_eps_max));
    cur_KL_start = max(1,sum(sorted_KL < KL_eps_min));
    
    sorted_KL = sorted_KL(cur_KL_start:cur_KL_end);
    Sanov_KL =Sanov_KL(cur_KL_start:cur_KL_end);


    % Now try to determine the scaling law according to Sanov theorem
    plot(sorted_KL, Sanov_KL+0*sorted_KL, colorvec(N_ind)); 
    title(['logprob/N KL 1-Cumulative Distribution for N = '  num2str(N_samples_vec(N_ind))]); xlabel('\epsilon'); ylabel('1/N log (1-F_{KL}(\epsilon))');
    
    
    
    
end % loop on N's

plot(sorted_KL, -sorted_KL, '--');

legend([num2str(N_samples_vec) ' lin-fit']);