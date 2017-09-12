cd output/linux;

% Exhaustive search
exh_score_linux = load('exhaust_learn_scores4_long.txt'); 
figure; hold on; plot( (res:res:max_nsamples), exh_scores ); title('Score of Exhaustive Searchs best DAG (N=4)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), exh_scores_linux, '.'); plot( (res:res:max_nsamples), exh_scores_linux, 'r'); 

% Max mcmc search
max_scores4_linux = load('mcmc_max_learn_scores4_long.txt');
figure; hold on; plot( (res:res:max_nsamples), max_score4(:,20) ); title('Score of MCMC Searchs best DAG (N=4)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), max_score4(:,1), '.'); plot( (res:res:max_nsamples), max_score4(:,1), 'r'); plot( (res:res:max_nsamples), max_scores4_linux(:,20) );
plot( (res:res:max_nsamples), max_scores4_linux(:,1), '.'); plot( (res:res:max_nsamples), max_scores4_linux(:,1), 'r'); legend('200 graphs', '10 graphs');
 
% Average mcmc search
ave_scores4_linux = load('mcmc_ave_learn_scores4_long.txt');
figure; hold on; plot( (res:res:max_nsamples), ave_score4(:,20) ); title('Score of MCMC Searchs average DAG (N=4)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), ave_score4(:,1), '.'); plot( (res:res:max_nsamples), ave_score4(:,1), 'r'); plot( (res:res:max_nsamples), ave_scores4_linux(:,20) );
plot( (res:res:max_nsamples), ave_scores4_linux(:,1), '.'); plot( (res:res:max_nsamples), ave_scores4_linux(:,1), 'r'); legend('200 graphs', '10 graphs');

% Weighted ave mcmc search
new_ave_scores4_linux = load('mcmc_new_ave_learn_scores4_long.txt');
figure; hold on; plot( (res:res:max_nsamples), new_ave_score4(:,20) ); title('Score of MCMC Searchs weighted ave DAG (N=4)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), new_ave_score4(:,1), '.'); plot( (res:res:max_nsamples), new_ave_score4(:,1), 'r'); plot( (res:res:max_nsamples), new_ave_scores4_linux(:,20) );
plot( (res:res:max_nsamples), new_ave_scores4_linux(:,1), '.'); plot( (res:res:max_nsamples), new_ave_scores4_linux(:,1), 'r'); legend('200 graphs', '10 graphs');



% Now Compare Methods ....
figure; hold on; plot( (res:res:max_nsamples), exh_scores ); plot( (res:res:max_nsamples), exh_scores_linux ); 
title('Comparison of all searches DAG (N=4)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), max_score4(:,20), 'r:' ); plot( (res:res:max_nsamples), max_scores4_linux(:,20), 'r:' ); 
plot( (res:res:max_nsamples), ave_score4(:,20), 'm-.' ); plot( (res:res:max_nsamples), ave_scores4_linux(:,20), 'm-.' ); 
plot( (res:res:max_nsamples), new_ave_score4(:,20), 'k--');  plot( (res:res:max_nsamples), new_ave_scores4_linux(:,20), 'k--'); 
legend('exhaust', 'max search', 'average search', 'weighted ave');



%%% Now the same for 5 nodes 
% Max mcmc search
max_scores5_linux = load('mcmc_max_learn_scores5_long.txt');
figure; hold on; plot( (res:res:max_nsamples), max_score5(:,20) ); title('Score of MCMC Searchs best DAG (N=5)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), max_score5(:,1), '.'); plot( (res:res:max_nsamples), max_score5(:,1), 'r'); plot( (res:res:max_nsamples), max_scores5_linux(:,20) );
plot( (res:res:max_nsamples), max_scores5_linux(:,1), '.'); plot( (res:res:max_nsamples), max_scores5_linux(:,1), 'r'); legend('200 graphs', '10 graphs');
 
% Average mcmc search
ave_scores5_linux = load('mcmc_ave_learn_scores5_long.txt');
figure; hold on; plot( (res:res:max_nsamples), ave_score5(:,20) ); title('Score of MCMC Searchs average DAG (N=5)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), ave_score5(:,1), '.'); plot( (res:res:max_nsamples), ave_score5(:,1), 'r'); plot( (res:res:max_nsamples), ave_scores5_linux(:,20) );
plot( (res:res:max_nsamples), ave_scores5_linux(:,1), '.'); plot( (res:res:max_nsamples), ave_scores5_linux(:,1), 'r'); legend('200 graphs', '10 graphs');

% Weighted ave mcmc search
new_ave_scores5_linux = load('mcmc_new_ave_learn_scores5_long.txt');
figure; hold on; plot( (res:res:max_nsamples), new_ave_score5(:,20) ); title('Score of MCMC Searchs weighted ave DAG (N=5)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), new_ave_score5(:,1), '.'); plot( (res:res:max_nsamples), new_ave_score5(:,1), 'r'); plot( (res:res:max_nsamples), new_ave_scores5_linux(:,20) );
plot( (res:res:max_nsamples), new_ave_scores5_linux(:,1), '.'); plot( (res:res:max_nsamples), new_ave_scores5_linux(:,1), 'r'); legend('200 graphs', '10 graphs');



% Now Compare Methods ....
figure; hold on; 
title('Comparison of all searches DAG (N=5)'); xlabel('num samples'); ylabel('errors score');  
plot( (res:res:max_nsamples), max_score5(:,20), 'r:' ); plot( (res:res:max_nsamples), max_scores5_linux(:,20), 'r:' ); 
plot( (res:res:max_nsamples), ave_score5(:,20), 'm-.' ); plot( (res:res:max_nsamples), ave_scores5_linux(:,20), 'm-.' ); 
plot( (res:res:max_nsamples), new_ave_score5(:,20), 'k--');  plot( (res:res:max_nsamples), new_ave_scores5_linux(:,20), 'k--'); 
legend('max search', 'average search', 'weighted ave');











%figure; surf(max_score4); title('4 nodes mc max score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
%figure; surf(ave_score4); title('4 nodes mc average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
%figure; surf(new_ave_score4); title('4 nodes mc new average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');



%figure; surf(max_score5); title('5 nodes mc max score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
%figure; surf(ave_score5); title('5 nodes mc average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');
%figure; surf(new_ave_score5); title('5 nodes mc new average score'); xlabel('graph samples'); ylabel('num samples'); zlabel('errors score');


cd ../..;