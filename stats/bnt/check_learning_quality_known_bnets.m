%%%function [max_learn_score, ave_learn_score, new_ave_score, max_exh_score, KL_max_learn_score, KL_ave_learn_score, KL_new_ave_score, KL_max_exh_score, correct_rank] = ...
function [edge_max_exh_score, KL_max_exh_score, edge_average_exh_score, correct_rank, edge_all_dags_exh_score, bayesian_all_dags_exh_score, KL_all_dags_exh_score, corr_bayesian_edge, corr_bayesian_KL] = ...
    check_learning_quality(points, nodes, max_nsamples, max_ngraph_samples, res, do_exh, do_mc, do_all_dags, do_KL, draw_fig)
% create 'random' bayesian num_points, and then try to learn them.
% We want to see how well we learn certain things 
% The number of input nodes is nodes.
% points are the probability functions defining the nets we are trying to
% learn


ep = 2.0/3.0;  % probability of an edge in the bnet    

directed = 1; undirected = 0;

num_points = size(points,2); % How many points


% Scores of edges
max_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
ave_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
new_ave_score = zeros(max_nsamples/res, max_ngraph_samples/res); 


% Scores of KL - relative entropy
KL_max_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
KL_ave_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
KL_new_ave_score = zeros(max_nsamples/res, max_ngraph_samples/res); 





edge_max_exh_score = zeros(num_points,max_nsamples/res);
KL_max_exh_score = zeros(num_points,max_nsamples/res);
edge_average_exh_score = zeros(num_points,max_nsamples/res);



% The place of the correct graph according to it's score
correct_rank = zeros(num_points, max_nsamples/res);



% generate the exhaustive DAGs
if(do_exh)
    exh_dags = zuk_mk_all_dags(nodes);
end

if(do_all_dags)
    edge_all_dags_exh_score  = zeros(num_points, length(exh_dags)); % Here we save the samples dimention ...
    bayesian_all_dags_exh_score = zeros(num_points, max_nsamples/res, length(exh_dags)); 
    KL_all_dags_exh_score  = zeros(num_points, max_nsamples/res, length(exh_dags));
end


ttt=cputime;
% Go over the number of bayesian nets needed

for i = 1:num_points
    
    
    %check time
    sprintf(' Finished Net %d\n', i)
    cputime-ttt
    ttt=cputime;
    
    % randomize the bayesian net
    %%%    bnet = create_random_bnet(nodes, ep);
    
    
    % Compute the bayesian net from the point
    ns = 2*ones(1,nodes); 
    
    dag = zeros(nodes);
    dag(1,3) = 1; dag(2,3) = 1;
    
    bnet = mk_bnet(dag, ns);
 
    CPT = [0.5 0.5];
    bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', CPT);
    bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', CPT);
    
    y = points(:,i);  
    CPT = [y(7),1-y(7); y(5),1-y(5);y(3),1-y(3);y(1),1-y(1)];
    bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', CPT);
    % Ended Compute the bayesian net from the point
    
    
    
    
    % Find the edge distances of all the dags
    if(do_all_dags)        
        for j=1:length(exh_dags)
            edge_all_dags_exh_score(i, j) =  dg_dag_dist(exh_dags{j}, bnet.dag, undirected);
        end
    end
        
    
 
    % Save the index of the correct dag
    for j=1:length(exh_dags)
        if(bnet.dag == exh_dags{j})
            correct_index = j;
            break;
        end
    end
    if(j == length(exh_dags))
       if( ~(bnet.dag == exh_dags{j}) )
           correct_index = -1;
       end
   end
   
    sprintf('The correct DAG is index %d\n', correct_index);
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Already inside the create_random_bnet    
    %derich_p = 1.0;
    
    % randomize the CPD's
    %for j=1:nodes
    %    ps = parents(bnet.dag, j);
    %    
    %    psz = power(2, length(ps)); %prod(ns(ps));
    %    CPT = dirichlet_sample(derich_p*ones(1,2), psz);
    %    bnet.CPD{j} = tabular_CPD(bnet, j, 'CPT', CPT);
    %
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    sam = cell(nodes, max_nsamples);
    
    % sample data from it
    for j = 1:max_nsamples
        sam(:,j) = sample_bnet(bnet);
    end
        
   % transfer from cell to a number matrix
   data = cell2num(sam);

    if(draw_fig)
         figure; draw_graph(bnet.dag); title('correct');
         figure; draw_graph(dag_to_cpdag(bnet.dag)); title('correct cpdag');
    end
   
   
   % Go and increase the samples, each time by res
   for(cur_samples = (res:res:max_nsamples))
        
       %sprintf('start i : %d\n', i)
       % Work only on part of the data
       cur_data = data(:, 1:cur_samples);     
   
       % Do exhaustive search
       if(nodes < 6)   
           if(do_exh)    
                %bayesian_exh_scores = score_dags(cur_data, 2*ones(1,nodes), exh_dags);
                bayesian_exh_scores = zuk_score_dags(cur_data, 2*ones(1,nodes), exh_dags);
%                exp(bayesian_exh_scores)
%                sprintf('SUM : %d \n', sum(exp(bayesian_exh_scores)))
                
                temp_bayesian_exh_scores = bayesian_exh_scores - max(bayesian_exh_scores);
                
                % Here calculate the 'average' DAG
                %temp_bayesian_exh_scores
                %sprintf('Now pause ..\n')
                %pause
                %exp(temp_bayesian_exh_scores)
                %sum(exp(temp_bayesian_exh_scores))
                bayesian_exh_probs = exp(temp_bayesian_exh_scores) / sum(exp(temp_bayesian_exh_scores));   
                %bayesian_exh_probs
                bayesian_average_dag = zeros(nodes);
                for k=1:length(exh_dags)
                    bayesian_average_dag = bayesian_average_dag + bayesian_exh_probs(k) * exh_dags{k};
                end
                bayesian_average_dag(find(bayesian_average_dag > 0.5)) = 1; bayesian_average_dag(find(bayesian_average_dag <= 0.5)) = 0;
                edge_average_exh_score(i, cur_samples/res) = dg_dag_dist(bayesian_average_dag, bnet.dag, undirected);
                
                
                % Find the index of the maximal bayesian score
                ex_max_index = find(bayesian_exh_scores == max(bayesian_exh_scores));
                edge_max_exh_score(i, cur_samples/res) = dg_dag_dist(exh_dags{ex_max_index(1)}, bnet.dag, undirected); % + edge_max_exh_score(i, cur_samples/res)
                
                if(draw_fig)
                    figure; draw_graph(exh_dags{ex_max_index(1)}); title('best exhaustive searched');
                end
                
                
                % Now find the 'place' of the correct dag
                correct_rank(i, cur_samples/res) = length( find(bayesian_exh_scores > bayesian_exh_scores(correct_index)) ) + 1;
                
                
                
                % Now learn parameters for the best bayesian scoring DAG
                temp_bnet = mk_bnet( exh_dags{ex_max_index(1)}, 2*ones(1,nodes));
                
                for k=1:nodes
                    temp_bnet.CPD{k} = tabular_CPD(temp_bnet, k);
                end
                
                temp_bnet2 = learn_params(temp_bnet, cur_data);    
                         
                KL_max_exh_score(i, cur_samples/res) =  zuk_bnet_relative_entropy(bnet, temp_bnet2); % + KL_max_exh_score(cur_samples/res) 
                
                
                
                
                
                if(do_all_dags)
                    bayesian_all_dags_exh_score(i, cur_samples/res, :) = bayesian_exh_scores;
                    
                    
                    if(do_KL)
                        % Enumerate all dags, learn parameters for each and
                        % score the DAGs according to this.
                        for j=1:length(exh_dags)
                            % Now learn parameters for the best bayesian scoring DAG
                            temp_bnet = mk_bnet( exh_dags{j}, 2*ones(1,nodes));
                    
                            for k=1:nodes
                                temp_bnet.CPD{k} = tabular_CPD(temp_bnet, k);
                            end
                
                            temp_bnet2 = learn_params(temp_bnet, cur_data);    
                         
                            KL_all_dags_exh_score(i, cur_samples/res, j) =  zuk_bnet_relative_entropy(bnet, temp_bnet2); % + KL_max_exh_score(cur_samples/res) 
                        end
                    end
                end
                
                
                
% %                 % Now learn parameters for each DAG
% %                 for i=1:length(exh_dags)
% %                      temp_bnet = mk_bnet(exh_dag(i), 2*ones(1:nodes));
% %                      temp_bnet2 = learn_params(temp_bnet, cur_data);
% % 
% %                      KL_exh_scores = zuk_bnet_relative_entropy(bnet, temp_bnet2);
% %                      
% %                      % Now calculate the Kullback-Leiber relative entropy 
% %                      KL_max_exh_score(cur_samples/res) = edge_max_exh_score(cur_samples/res)  + 
% %     
% %                 end    
% %                 
% %                 KL_ex_max_index = find(KL_exh_scores == max(KL_exh_scores));
% %                 KL_max_exh_score(cur_samples/res) = KL_max_exh_score(cur_samples/res) + dg_dag_dist(exh_dags{ex_max_index(1)}, bnet.dag, undirected);
                
           end
       end 
           
    end % loop over nsamples
    
%    sprintf('Learned Score : %d\n', sum(sum(qu{max_index(1)} == bnet.dag)))
%    sprintf('AveLear Score : %d\n', sum(sum(ave_dag == bnet.dag)))
   
   
end % loop over nets 



corr_bayesian_edge = zeros(num_points, max_nsamples/res);
corr_bayesian_KL = zeros(num_points, max_nsamples/res);
% Calculate correlations
if(do_all_dags)
    for i = 1:num_points
        for j = 1:max_nsamples/res
%             sprintf('The Sizes (edge : bayes) \n')
%             size(edge_all_dags_exh_score)
%             size(bayesian_all_dags_exh_score)
            
            temp = corrcoef(edge_all_dags_exh_score(i,:), bayesian_all_dags_exh_score(i, j, :));
            corr_bayesian_edge(i, j) = temp(2);
        
            if(do_KL)
%                 sprintf('The Sizes (KL : bayes)\n')
%                 size(KL_all_dags_exh_score)
%                 size(bayesian_all_dags_exh_score)
                temp = corrcoef(KL_all_dags_exh_score(i, j, :), bayesian_all_dags_exh_score(i, j, :));
                corr_bayesian_KL(i, j) = temp(2);
            end
        end
    end

% % %     edge_all_dags_exh_score  = zeros(num_points, length(exh_dags)); % Here we save the samples dimention ...
% % %     bayesian_all_dags_exh_score = zeros(num_points, max_nsamples/res, length(exh_dags)); 
% % %     KL_all_dags_exh_score  = zeros(num_points, max_nsamples/res, length(exh_dags));
end    
    



   
%bnet.dag;
%dag_to_cpdag(bnet.dag);
%qu{max_index(1)};
%new_ave_dag;
%exh_dags{ex_max_index(1)};

%sprintf('Learned Score : %d\n', max_learn_score);
%sprintf('AveLear Score : %d\n', ave_learn_score);
%sprintf('NewAve  Score : %d\n', new_ave_score);
if(nodes < 6)
    if(do_exh)
        sprintf('Exhaustive Search Score : %d\n', edge_max_exh_score);
    end
end
%gr = sgu; %sampled_graphs;


% Here we return things 
max_learn_score = max_learn_score/num_points; 
ave_learn_score = ave_learn_score/num_points;
new_ave_score = new_ave_score/num_points;
%edge_max_exh_score = edge_max_exh_score/num_points;

KL_max_learn_score = KL_max_learn_score/num_points; 
KL_ave_learn_score = KL_ave_learn_score/num_points;
KL_new_ave_score = KL_new_ave_score/num_points;
%KL_max_exh_score = KL_max_exh_score/num_points;



plot(correct_rank);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare two DGs to see how many edges differ. The first is the directed
% graph we've gor. The second is a DAG - usually the original DAG that we
% sampled from. If directed is on, we account also for the direction of
% edges learned. Otherwise we account only for edges/non edges
function dd = dg_dag_dist(dg1, dag2, directed)

dg1;
dag2;
%dag_to_cpdag(dag1)
%dag_to_cpdag(dag2)

% first transfer into pdags
%pdag1 = dag_to_cpdag(dag1) + dag1;  % 1 for undirected edge (in both places), 2 for directed edge (in one place ..)
%pdag2 = dag_to_cpdag(dag2) + dag2;

nodes = size(dag2, 1);

% Here the undirected case ...
if(directed == 0)
     un_dg1 = dg1+dg1';
     un_dag2 = dag2+dag2';
%     un_dg1 .* un_dag2
%     find(un_dg1 .* un_dag2)
%     length(find(un_dg1 .* un_dag2))
    dd = ( length(find(un_dg1 + un_dag2)) - length(find(un_dg1 .* un_dag2)) ) / (nodes * (nodes-1));
else

    %sprintf('the two compressed :\n')
    comp_dg1 = digraph_to_compressed_digraph(dg1);
    comp_dag2 = digraph_to_compressed_digraph(dag_to_cpdag(dag2));




    % 0 - no edge at all
    % 1 - edge from i to j
    % 2 - edge from j to i
    % 3 - undirected edge i-j (both directions are present).
    % Now calc the distance

    % Now calculate the score . Punishments are : 
    % 0 '   '  0  '   '    0
    % 3 '---'  3  '---'    0
    % 1 '-->'  1  '-->'    0
    % 0 '   '  3 '---'     3
    % 0 '   '  1 '-->'     3
    % 3 '---'  1 '-->'     1
    % 1 '-->'  2 '<--'     2

                  % 0  1  2  3
    new_pun_tab = [ 0, 3, 3, 3; ... %0
                    3, 0, 2, 1; ... %1
                    3, 2, 0, 1; ... %2
                    3, 1, 1, 0];    %3 

        
    % Non symetrical punishment table. 
    % This is because in the learned DAG, we cannot
    % excpect to find both directions of an undirected edge
    % Note : The i is for the learned, the j for the original
                       % 0  1  2  3
    new_pun_tab_asym = [ 0, 3, 3, 3; ... %0
                         3, 0, 2, 0; ... %1
                         3, 2, 0, 0; ... %2
                         3, 1, 1, 0];    %3        

    punish = 0;


    for i = (1:nodes)
        for j = (i+1:nodes)
         % sprintf('dg1 : %d dag2 : %d\n', comp_dg1(i, j), comp_dag2(i, j))
            punish = punish + new_pun_tab_asym(comp_dg1(i, j)+1, comp_dag2(i, j)+1);    
        end
    end

    % We have to divide it by something, so to denote the 'precentage' of wrong
    % edges
    dd = punish/ ( 3*nodes*(nodes-1)/2 );
        
end % Of else part


        






























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   
% just to save ....
function dd = dg_dag_dist_old(dg1, dag2)

cpdag1 = dag_to_cpdag(dag1);
pdag1 = cpdag1 .* ( 3 - (cpdag1 + cpdag1') );
cpdag2 = dag_to_cpdag(dag2);
pdag2 = cpdag2 .* ( 3 - (cpdag2 + cpdag2') );


pdag1;
pdag2;


nodes = size(dag1);


% prepare the punish tab
%pun_tab = zeros(3, 3, 3, 3);


          % 00 02 20 11 
pun_tab = [ 0, 3, 3, 3; ... % 00 
            3, 0, 2, 1; ... % 02
            3, 2, 0, 1; ... % 20 
            3, 1, 1, 0];   % 11
        
        
        
                



% Now calculate the score . Punishments are : 
% '   '    '   '    0
% '---'    '---'    0
% '-->'    '-->'    0
% '   '    '---'    3
% '   '    '-->'    3
% '---'    '-->'    1
% '-->'    '<--'    2
punish = 0;

for i = [1:nodes]
    for j = [i+1:nodes]
        sprintf('IJ %d %d %d %d ind %d %d', pdag1(i, j), pdag1(j, i), pdag2(i, j), pdag2(j, i), ...
                                  1 +  pdag1(i,j) + pdag1(j, i)/2 + 3*pdag1(i,j)*pdag1(j, i)/2,  ...
                                  1 +  pdag2(i,j) + pdag2(j, i)/2 + 3*pdag2(i,j)*pdag2(j, i)/2)
        punish = punish + pun_tab(1 +  pdag1(i,j) + pdag1(j, i)/2 + 3*pdag1(i,j)*pdag1(j, i)/2,  ...
                                  1 +  pdag2(i,j) + pdag2(j, i)/2 + 3*pdag2(i,j)*pdag2(j, i)/2);
    end
end

dd = punish;

