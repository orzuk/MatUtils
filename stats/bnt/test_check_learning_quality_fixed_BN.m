% determine parameters
% path(path, 'E:\Research\networks\BNsoftware\KevinMurphy\FullBNT\BNT\general');
% path(path, 'E:\Research\networks\BNsoftware\KevinMurphy\FullBNT\graph');
% path(path, 'E:\Research\networks\BNsoftware\KevinMurphy\FullBNT\KPMtools');
% path(path, 'E:\Research\networks\BNsoftware\KevinMurphy\FullBNT\KPMstats');

nets = 20;  % How many times to randomize a BNT
nodes = 4; 

iters = 100; % How many iterations for each network
max_nsamples = 500;
res = 100;
%ngraph_samples = 100;
max_ngraph_samples = res;

epsilon = 0.00000001;

draw_fig = 0; % 0
do_exh = 1; do_mc = 0; % Exhaustive search all DAGs or sample DAGs.
do_save = 0; do_log = 0;  % Save to file or not
do_all_dags = 1;   % 1 - score all DAGs exhaustively. 0 - Score only best DAG.
do_run = 1;    % 1 - run everything. 0 - only plot graphs.
do_KL = 0;     % 1 - compute KL scores for all the DAGs (only if do_all is on). 0 - don't compute.


ep = 0.5; % Probability for something

if(do_run)
    sprintf('Start Exh 4 ')
    % start learning 
    %%%edge_exh_scores = zeros(1, max_nsamples/res); KL_exh_scores = zeros(1, max_nsamples/res); 
    ttt=cputime;
    

    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Generate my BNET 
     %%%bnet = create_random_bnet(nodes, ep);  % How do we randomize ????? 
    
     
     
     % create random DAG
    dag = zeros(nodes); dag(1,2)=1; dag(1,3)=1; dag(2,4)=1; dag(3,4)=1;
    
%%%    dag = zeros(nodes); dag(1,2)=1; dag(1,3)=1; %dag(2,4)=1; 
    
    
   

% Now make a bayesian networks from the dag. 
% For now, we assume all nodes are discrete and boolean
ns = 2*ones(1,nodes); 

bnet = mk_bnet(dag, ns);


% % Now account for the CPD's
    CPT = [0.9 0.1];
    bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', CPT);

    CPT = [0.9 0.1; 0.7 0.3];
     bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', CPT);
    
    CPT = [0.9 0.1; 0.7 0.3];
    bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', CPT);
    

    CPT = [0.9 0.1; 0.7 0.3; 0.2 0.8; 0.8 0.2];
    bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', CPT);
%     
% % Finished generating my BNET
% %     CPT = [0.6 0.5];
% %     bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', CPT);
% % 
% %     CPT = [0.55 0.45; 0.45 0.55];
% %      bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', CPT);
% %     
% %     CPT = [0.2 0.8; 0.43 0.57];
% %     bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', CPT);
% %     
% % 
% %     CPT = [0.35 0.65]; %%%  0.27 0.73];
% %     bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', CPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
    
    
    %[dum_score, dum_ave_score, dum_new_ave_score, exh_scores, KL_dum_score, KL_dum_ave_score, KL_dum_new_ave_score, KL_exh_scores] = ...
    [ edge_exh_scores4, KL_exh_scores4, edge_average_exh_scores4, correct_ranks4, edge_all_dags_exh_scores4, ...
        bayesian_all_dags_exh_scores4, KL_all_dags_exh_scores4, corr_bayesian_edge4, corr_bayesian_KL4, exh_dags, ...
         dags_order_mat, containing_error, not_containing_error] = ...
        check_learning_quality_fixed_BN(bnet, nodes, iters,  max_nsamples, max_ngraph_samples, res, do_exh, do_mc, do_all_dags, do_KL, draw_fig);
    alll_time_passed = cputime-ttt
end % do_run



if(do_save)
    dir_name = strcat('nodes.', num2str(nodes), 'nets.', num2str(nets), 'samples.', num2str(max_nsamples));
    mkdir(sprintf('%s', dir_name));
    cd(sprintf('%s', dir_name));
    
    
    % Save to a file    
    save 'corr_bayesian_edge.txt' corr_bayesian_edge4 -ASCII;
    save 'exhaust_diff_edge_scores.txt' edge_exh_scores4 -ASCII;
    save 'exhaust_diff_KL_scores.txt' KL_exh_scores4 -ASCII;
    save 'exhaust_diff_edge_average_scores.txt' edge_average_exh_scores4 -ASCII;
    
    if(do_KL)
        save 'corr_bayesian_KL.txt' corr_bayesian_edge4 -ASCII;
    end 
end






return; % Avoid all annoying figures









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Start Plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do a readable subplot of three for Shiri's thesis & paper etc.
if(do_log)
    KL_exh_scores4 = log(KL_exh_scores4 + epsilon);
end
mean_KL_exh_scores4 = mean(KL_exh_scores4, 1);
std_KL_exh_scores4 = std(KL_exh_scores4, 1);
h = figure; hold on; subplot(3,1,1); hold on;  % The KL distance
title('KL Dist of correct and learned DAG'); xlabel('Samples Num. '); ylabel('KL Dist.'); 
errorbar((res:res:max_nsamples), mean_KL_exh_scores4, std_KL_exh_scores4, std_KL_exh_scores4);


subplot(3,1,2); hold on; % Rank of correct struct.
title('Rank of the correct structure among all DAGs'); xlabel('Samples Num. '); ylabel('Rank'); 
mean_correct_ranks4 = mean(correct_ranks4, 1);
std_correct_ranks4 = std(correct_ranks4, 1);
errorbar((res:res:max_nsamples), mean_correct_ranks4, std_correct_ranks4, std_correct_ranks4);

subplot(3,1,3); hold on; % Fraction of correctly recovered edges.
title('Edge Errors of learned DAG'); xlabel('Samples Num. '); ylabel('Diff Edges Fraction');  
if(do_log)
    edge_exh_scores4 = log(edge_exh_scores4 + epsilon);
end
mean_edge_exh_scores4 = mean(edge_exh_scores4, 1);
std_edge_exh_scores4 = std(edge_exh_scores4, 1);
errorbar((res:res:max_nsamples), mean_edge_exh_scores4, std_edge_exh_scores4, std_edge_exh_scores4); 









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Seperate figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the edge score of best DAG as function of samples
h = figure; hold on; 
title('Edge Score of Exhaustive Searchs best DAG 4'); xlabel('num samples'); ylabel('errors score');  
if(do_log)
    edge_exh_scores4 = log(edge_exh_scores4 + epsilon);
end
mean_edge_exh_scores4 = mean(edge_exh_scores4, 1);
std_edge_exh_scores4 = std(edge_exh_scores4, 1);
errorbar((res:res:max_nsamples), mean_edge_exh_scores4, std_edge_exh_scores4, std_edge_exh_scores4); 

% % % if(do_save)
% % %     saveas(h, 'exhaust_diff_edge_scores.fig');
% % % end


% New ! Plot everything in log-log plot
h = figure; 
if(do_log)
    KL_exh_scores4 = log(KL_exh_scores4 + epsilon);
end
mean_KL_exh_scores4 = mean(KL_exh_scores4, 1);
std_KL_exh_scores4 = std(KL_exh_scores4, 1);
loglog(res:res:max_nsamples, mean_KL_exh_scores4); hold on; 
title('LogLog KL Score of Exhaustive Searchs best DAG 4'); xlabel('num samples'); ylabel('errors score'); 
% Plot the KL score of best DAG as function of samples
h = figure; hold on;
title('KL Score of Exhaustive Searchs best DAG 4'); xlabel('num samples'); ylabel('errors score'); 
errorbar((res:res:max_nsamples), mean_KL_exh_scores4, std_KL_exh_scores4, std_KL_exh_scores4);






% % % if(do_save)
% % %     saveas(h, 'exhaust_diff_KL_scores.fig');
% % % end


% Plot the edge score of average DAG from bayesian averaging as function of samples
h = figure; hold on; 
title('Edge Score of Exhaustive Searchs average DAG 4'); xlabel('num samples'); ylabel('errors score');  
if(do_log)
    edge_average_exh_scores4 = log(edge_average_exh_scores4 + epsilon);
end
mean_edge_average_exh_scores4 = mean(edge_average_exh_scores4, 1);
std_edge_average_exh_scores4 = std(edge_average_exh_scores4, 1);
errorbar((res:res:max_nsamples), mean_edge_average_exh_scores4, std_edge_average_exh_scores4, std_edge_average_exh_scores4); 



% Plot the rank of the correct DAG according to the bayesian score as
% function of number of samples
h = figure; hold on;
title('Rank of the correct structure of Exhaustive Searchs best DAG 4'); xlabel('num samples'); ylabel('rank'); 
mean_correct_ranks4 = mean(correct_ranks4, 1);
std_correct_ranks4 = std(correct_ranks4, 1);
errorbar((res:res:max_nsamples), mean_correct_ranks4, std_correct_ranks4, std_correct_ranks4);
if(do_save)
    
    % Save to a file    
    save 'correct_ranks.txt' edge_exh_scores4 -ASCII;    
    % % %     saveas(h, 'correct_ranks.fig');
end


% Plot correlation between the the edge score and KL score of the best DAG as a function of the number of samples 
h = figure; hold on;
title('Correlation between edge score and KL score of Exhaustive Searchs best DAG 4'); xlabel('num samples'); ylabel('correlation'); 
corr4 = zeros(1, max_nsamples/res);
for j=1:(max_nsamples/res)
    temp = corrcoef(edge_exh_scores4(:,j), KL_exh_scores4(:,j));  
    corr4(j) = temp(2);
end
length(corr4)
length(res:res:max_nsamples)
plot((res:res:max_nsamples), corr4, '.');

if(do_save)
    save 'corr_edge_KL_exhaust.txt' corr4 -ASCII;
end

% Now do some global plots - not only the best scoring one
if(do_all_dags)
    want_to_do_many_plots = 0; 
    if(want_to_do_many_plots)
    
    % Plot correlation between edge score and bayesian score
    mean_corr_bayesian_edge4 = mean(corr_bayesian_edge4, 1);
    std_corr_bayesian_edge4 = std(corr_bayesian_edge4, 1);
    figure; hold on; errorbar((res:res:max_nsamples), mean_corr_bayesian_edge4, std_corr_bayesian_edge4, std_corr_bayesian_edge4);
    xlabel('num samples'); ylabel('correlation'); title('correlation between bayesian score and edge score of all dags, DAG 4');
    
    
    
    % Do some sample plottings of all the DAGS ..
    i=1;
    figure; hold on; plot(reshape( bayesian_all_dags_exh_scores4(1, i, :), [1 446] ), edge_all_dags_exh_scores4(1, :), '.'); 
    xlabel('bayesian score'); ylabel('edge score'); title(strcat('one net plot of correctness (edge) vs. bayesian score  ', num2str(i*res), ' samples'));
    
    i=max_nsamples/res;
    figure; hold on; plot(reshape( bayesian_all_dags_exh_scores4(1, i, :), [1 446] ), edge_all_dags_exh_scores4(1, :), '.'); 
    xlabel('bayesian score'); ylabel('edge score'); title(strcat('one net plot of correctness (edge) vs. bayesian score  ', num2str(i*res), ' samples'));
    
    % No do averaging on all the nets
    figure; hold on; 
    xlabel('bayesian score'); ylabel('edge score'); title(strcat('all nets plot of correctness (edge) vs. bayesian score  ', num2str(i*res), ' samples'));
    for i = 1:max_nsamples/res
        plot(reshape( bayesian_all_dags_exh_scores4(1, i, :), [1 446] ), edge_all_dags_exh_scores4(1, :), '.'); 
    end 
    
    
    
    if(do_KL)
        % Plot correlation between KL score and bayesian score
        mean_corr_bayesian_KL4 = mean(corr_bayesian_KL4, 1);
        std_corr_bayesian_KL4 = std(corr_bayesian_KL4, 1);
        figure; hold on; errorbar((res:res:max_nsamples), mean_corr_bayesian_KL4, std_corr_bayesian_KL4, std_corr_bayesian_KL4);
        xlabel('num samples'); ylabel('correlation'); title('correlation between KL score and edge score of all dags, DAG 4');       
        
    end    
    
end
    
end    








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















cd ..
    sprintf('Before Start Exh 5 ')
while( 1==0) % Never start  ???
    
    % Do the same for 5 nodes - might be very slow !!!!!!!!!
    nodes = 5; 
    
    sprintf('Start Exh 5 ')
    % start learning 
    %%%edge_exh_scores = zeros(1, max_nsamples/res); KL_exh_scores = zeros(1, max_nsamples/res); 
    ttt=cputime;
    %[dum_score, dum_ave_score, dum_new_ave_score, edge_exh_scores, KL_dum_score, KL_dum_ave_score, KL_dum_new_ave_score, KL_exh_scores] = ...
    [ edge_exh_scores5, KL_exh_scores5, edge_average_exh_scores5, correct_ranks5, edge_all_dags_exh_scores5, bayesian_all_dags_exh_scores5, KL_all_dags_exh_scores5,  corr_bayesian_edge5, corr_bayesian_KL5] = ...
        check_learning_quality(nets, nodes, max_nsamples, max_ngraph_samples, res, do_exh, do_mc, do_all_dags, draw_fig);
    cputime-ttt
    
    if(do_save)
        dir_name = strcat('nodes.', num2str(nodes), 'nets.', num2str(nets), 'samples.', num2str(max_nsamples));
        cd ..;
        mkdir(sprintf('%s', dir_name));
        cd(sprintf('%s', dir_name));
        
        % Save to a file    
        save 'corr_bayesian_edge.txt' corr_bayesian_edge5 -ASCII;
        save 'exhaust_diff_edge_scores.txt' edge_exh_scores5 -ASCII;
        save 'exhaust_diff_KL_scores.txt' KL_exh_scores5 -ASCII;
        
        if(do_KL)
            save 'corr_bayesian_KL.txt' corr_bayesian_edge5 -ASCII;
        end 
        
    end
    
    
    
    % Plot the edge score of best DAG as function of samples
    h = figure; hold on; 
    title('Edge Score of Exhaustive Searchs best DAG 5'); xlabel('num samples'); ylabel('errors score');  
    if(do_log)
        edge_exh_scores5 = log(edge_exh_scores5 + epsilon);
    end
    mean_edge_exh_scores5 = mean(edge_exh_scores5, 1);
    std_edge_exh_scores5 = std(edge_exh_scores5, 1);
    errorbar((res:res:max_nsamples), mean_edge_exh_scores5, std_edge_exh_scores5, std_edge_exh_scores5); 
    
    % % % if(do_save)
    % % %     saveas(h, 'exhaust_diff_edge_scores.fig');
    % % % end
    
    
    
    % Plot the KL score of best DAG as function of samples
    h = figure; hold on;
    title('KL Score of Exhaustive Searchs best DAG 5'); xlabel('num samples'); ylabel('errors score'); 
    if(do_log)
        KL_exh_scores5 = log(KL_exh_scores5 + epsilon);
    end
    mean_KL_exh_scores5 = mean(KL_exh_scores5, 1);
    std_KL_exh_scores5 = std(KL_exh_scores5, 1);
    errorbar((res:res:max_nsamples), mean_KL_exh_scores5, std_KL_exh_scores5, std_KL_exh_scores5);
    
    
    
    
    % % % if(do_save)
    % % %     saveas(h, 'exhaust_diff_KL_scores.fig');
    % % % end
    
    
    % Plot the rank of the correct DAG according to the bayesian score as
    % function of number of samples
    h = figure; hold on;
    title('Rank of the correct structure of Exhaustive Searchs best DAG 5'); xlabel('num samples'); ylabel('rank'); 
    mean_correct_ranks5 = mean(correct_ranks5, 1);
    std_correct_ranks5 = std(correct_ranks5, 1);
    errorbar((res:res:max_nsamples), mean_correct_ranks5, std_correct_ranks5, std_correct_ranks5);
    if(do_save)
        
        % Save to a file    
        save 'correct_ranks.txt' edge_exh_scores5 -ASCII;    
        % % %     saveas(h, 'correct_ranks.fig');
    end
    
    
    % Plot correlation between the the edge score and KL score of the best DAG as a function of the nuber of samples 
    h = figure; hold on;
    title('Corrlation between edge score and KL score of Exhaustive Searchs best DAG 5'); xlabel('num samples'); ylabel('correlation'); 
    corr5 = zeros(1, max_nsamples/res);
    for i=1:(max_nsamples/res)
        temp = corrcoef(edge_exh_scores5(:,i), KL_exh_scores5(:,i));  
        corr5(i) = temp(2);
    end
    length(corr5)
    length(res:res:max_nsamples)
    plot((res:res:max_nsamples), corr5);
    
    if(do_save)
        save 'corr_edge_KL_exhaust.txt' corr5 -ASCII;
    end
    
    
    
    
    
    
    
    
    % plot 
    % % % figure; plot( (res:res:max_nsamples), edge_exh_scores5 ); title('Edge Score of Exhaustive Searchs best DAG 5'); xlabel('num samples'); ylabel('errors score');  
    % % % figure; plot( (res:res:max_nsamples), KL_exh_scores5, 'r' );  title('KL Score of Exhaustive Searchs best DAG 5'); xlabel('num samples'); ylabel('errors score');  
    
    
    
    
    %while(1==0)
    
    do_exh = 0; do_mc = 1;
    
    
    sprintf('Start MCMC 4  ')
    
    max_score4 = zeros(max_nsamples, max_ngraph_samples);
    ave_score4 = zeros(max_nsamples, max_ngraph_samples);
    new_ave_score4 = zeros(max_nsamples, max_ngraph_samples);
    KL_max_score4 = zeros(max_nsamples, max_ngraph_samples);
    KL_ave_score4 = zeros(max_nsamples, max_ngraph_samples);
    KL_new_ave_score4 = zeros(max_nsamples, max_ngraph_samples);
    
    % Now learn non-exhaustively using the mcmc 
    %for nsamples =(res:res:max_nsamples)
    %    for ngraph_samples = (res:res:max_ngraph_samples)
    [max_score4, ave_score4, new_ave_score4, dum_exh_score, KL_max_score4, KL_ave_score4, KL_new_ave_score4, KL_dum_exh_score] = ...
        check_learning_quality(nets, nodes, max_nsamples, max_ngraph_samples, res, do_exh, do_mc, draw_fig);
    %    end
    %end
    
    % Save to a file
    save 'mcmc_max_learn_scores4_new.txt' max_score4 -ASCII;
    save 'mcmc_ave_learn_scores4_new.txt' ave_score4 -ASCII;
    save 'mcmc_new_ave_learn_scores4_new.txt' new_ave_score4 -ASCII;
    
    
    figure; surf(max_score4); title('4 nodes mc max score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
    figure; surf(ave_score4); title('4 nodes mc average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
    figure; surf(new_ave_score4); title('4 nodes mc new average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
    
    
    
    sprintf('Start MCMC 5  ')
    
    
    % Now do the same with five nodes
    nodes = 5; 
    max_score5 = zeros(max_nsamples, max_ngraph_samples);
    ave_score5 = zeros(max_nsamples, max_ngraph_samples);
    new_ave_score5 = zeros(max_nsamples, max_ngraph_samples);
    KL_max_score5 = zeros(max_nsamples, max_ngraph_samples);
    KL_ave_score5 = zeros(max_nsamples, max_ngraph_samples);
    KL_new_ave_score5 = zeros(max_nsamples, max_ngraph_samples);
    
    % Now learn non-exhaustively using the mcmc 
    % for nsamples =(res:res:max_nsamples)
    %     for ngraph_samples = (res:res:max_ngraph_samples)
    [max_score5, ave_score5, new_ave_score5, dum_exh_score, KL_max_score5, KL_ave_score5, KL_new_ave_score5, KL_dum_exh_score] = ...
        check_learning_quality(nets, nodes, max_nsamples, max_ngraph_samples, res, do_exh, do_mc, draw_fig);
    %     end
    % end
    
    
    % Save to a file
    save 'mcmc_max_learn_scores5_new.txt' max_score5 -ASCII;
    save 'mcmc_ave_learn_scores5_new.txt' ave_score5 -ASCII;
    save 'mcmc_new_ave_learn_scores5_new.txt' new_ave_score5 -ASCII;
    
    
    figure; surf(max_score5); title('5 nodes mc max score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
    figure; surf(ave_score5); title('5 nodes mc average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
    figure; surf(new_ave_score5); title('5 nodes mc new average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
    
    
end % artificial 'while' 