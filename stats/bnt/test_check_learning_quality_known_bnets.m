% determine parameters

path(path, 'C:\Weizmann\Research\BayesianNetworks\BNsoftware\KevinMurphy\FullBNT_2005_01_31\FullBNT\graph');
path(path, 'C:\Weizmann\Research\BayesianNetworks\BNsoftware\KevinMurphy\FullBNT_2005_01_31\FullBNT\BNT\CPDs\@tabular_CPD');

nodes = 3; 

max_nsamples = 500;
res = 500;
%ngraph_samples = 100;
max_ngraph_samples = 5;

epsilon = 0.00000001;

draw_fig = 0;
do_exh = 1; do_mc = 0;
do_save = 0; do_log = 0;  % Save to file or not
do_all_dags = 1;   % 1 - score all DAGs exhaustively. 0 - Score only best DAG.
do_run = 1;    % 1 - run everything. 0 - only plot graphs.
do_KL = 0;     % 1 - compute KL scores for all the DAGs (only if do_all is on). 0 - don't compute.


if(do_run)
sprintf('Start Exh 4 ')
% start learning 
%%%edge_exh_scores = zeros(1, max_nsamples/res); KL_exh_scores = zeros(1, max_nsamples/res); 
points_eps = 1;
num_points = 20;
ttt=cputime;
    points = generate_points(num_points,points_eps)


%[dum_score, dum_ave_score, dum_new_ave_score, exh_scores, KL_dum_score, KL_dum_ave_score, KL_dum_new_ave_score, KL_exh_scores] = ...
    [ edge_exh_scores3, KL_exh_scores3, edge_average_exh_scores3, correct_ranks3, edge_all_dags_exh_scores3, bayesian_all_dags_exh_scores3, KL_all_dags_exh_scores3, corr_bayesian_edge3, corr_bayesian_KL3] = ...
    check_learning_quality_known_bnets(points, nodes, max_nsamples, max_ngraph_samples, res, do_exh, do_mc, do_all_dags, do_KL, draw_fig);
cputime-ttt
end % do_run





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Start Plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the edge score of best DAG as function of samples
% % % h = figure; hold on; 
% % % title('Edge Score of Exhaustive Searchs best DAG 3'); xlabel('num samples'); ylabel('errors score');  
% % % if(do_log)
% % %     edge_exh_scores3 = log(edge_exh_scores3 + epsilon);
% % % end
% % % mean_edge_exh_scores3 = mean(edge_exh_scores3, 1);
% % % std_edge_exh_scores3 = std(edge_exh_scores3, 1);
% % % errorbar((res:res:max_nsamples), mean_edge_exh_scores3, std_edge_exh_scores3, std_edge_exh_scores3); 
% % % 
% % % % % % if(do_save)
% % % % % %     saveas(h, 'exhaust_diff_edge_scores.fig');
% % % % % % end
% % % 
% % % 
% % % 
% % % % Plot the KL score of best DAG as function of samples
% % % h = figure; hold on;
% % % title('KL Score of Exhaustive Searchs best DAG 3'); xlabel('num samples'); ylabel('errors score'); 
% % % if(do_log)
% % %     KL_exh_scores3 = log(KL_exh_scores3 + epsilon);
% % % end
% % % mean_KL_exh_scores3 = mean(KL_exh_scores3, 1);
% % % std_KL_exh_scores3 = std(KL_exh_scores3, 1);
% % % errorbar((res:res:max_nsamples), mean_KL_exh_scores3, std_KL_exh_scores3, std_KL_exh_scores3);
% % % 
% % % 
% % % % % % if(do_save)
% % % % % %     saveas(h, 'exhaust_diff_KL_scores.fig');
% % % % % % end
% % % 
% % % 
% % % % Plot the edge score of average DAG from bayesian averaging as function of samples
% % % h = figure; hold on; 
% % % title('Edge Score of Exhaustive Searchs average DAG 3'); xlabel('num samples'); ylabel('errors score');  
% % % if(do_log)
% % %     edge_average_exh_scores3 = log(edge_average_exh_scores3 + epsilon);
% % % end
% % % mean_edge_average_exh_scores3 = mean(edge_average_exh_scores3, 1);
% % % std_edge_average_exh_scores3 = std(edge_average_exh_scores3, 1);
% % % errorbar((res:res:max_nsamples), mean_edge_average_exh_scores3, std_edge_average_exh_scores3, std_edge_average_exh_scores3); 



% Plot the rank of the correct DAG according to the bayesian score as
% function of number of samples
h = figure; hold on;
title('Rank of the correct structure for different bnets (epsilons)'); xlabel('epsilon'); ylabel('rank'); 
plot(correct_ranks3,'.');


cd C:\Weizmann\Research\BayesianNetworks\BNsoftware\KevinMurphy\FullBNT_2005_01_31\FullBNT\shiri;

