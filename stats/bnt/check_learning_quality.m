%%%function [max_learn_score, ave_learn_score, new_ave_score, max_exh_score, KL_max_learn_score, KL_ave_learn_score, KL_new_ave_score, KL_max_exh_score, correct_rank] = ...
%%% New ! Changes made from 31.2.05 (Zuk)
function [edge_max_exh_score, KL_max_exh_score, edge_average_exh_score, correct_rank, edge_all_dags_exh_score, bayesian_all_dags_exh_score, KL_all_dags_exh_score, corr_bayesian_edge, corr_bayesian_KL] = ...
    check_learning_quality(nets, nodes, max_nsamples, max_ngraph_samples, res, do_exh, do_mc, do_all_dags, do_KL, draw_fig)
% create 'random' bayesian nets, and then try to learn them.
% We want to see how well we learn certain things 
% The number of input nodes is nodes.
% nets is the number of iterations - how many nets we try to learn.
% max_nsamples is the maximal number of samples we draw for each graph
% max_ngraph_samples is i don't have the slightest idea what is it ????
% res is the resulution in which we increment the of samples each time.
% do_exh is a flag saying if we enumerate all graphs exhaustively or just
% sample the space of DAGs
% do_mc is a flag saying to do monte-carlo search on the space of DAGs
% (???)
% do_all_dags is a flag saying if we do all the DAGs  - doh
% do_KL is a flag saying if we want to compute the KL distance between the
% correct and learned graph.
% draw_fig is a flag saying if we want to draw zillion figures or not



ep = 2.0/3.0;  % probability of an edge in the bnet    

directed = 1; undirected = 0; 
only_connected = 0;  % flag saying if we want only connected graphs or all graphs ! 


% Scores of edges
max_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
ave_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
new_ave_score = zeros(max_nsamples/res, max_ngraph_samples/res); 

% Scores of KL - relative entropy
KL_max_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
KL_ave_learn_score = zeros(max_nsamples/res, max_ngraph_samples/res); 
KL_new_ave_score = zeros(max_nsamples/res, max_ngraph_samples/res); 


edge_max_exh_score = zeros(nets,max_nsamples/res);
KL_max_exh_score = zeros(nets,max_nsamples/res);
edge_average_exh_score = zeros(nets,max_nsamples/res);


% The place of the correct graph according to it's score
correct_rank = zeros(nets, max_nsamples/res);



% generate the exhaustive DAGs
if(do_exh)
    exh_dags = zuk_mk_all_dags(nodes, only_connected);
end



% New : this vector says if the true probability distribution is contained
% in the DAG 
iss_contained_flags = zeros(1, length(exh_dags));



% Here loop over all Dags and ... Default for do_all_dags is one !
if(do_all_dags)
    edge_all_dags_exh_score  = zeros(nets, length(exh_dags)); % Here we save the samples dimention ...
    bayesian_all_dags_exh_score = zeros(nets, max_nsamples/res, length(exh_dags)); 
    KL_all_dags_exh_score  = zeros(nets, max_nsamples/res, length(exh_dags));
end


ttt=cputime;



% Go over the number of bayesian nets needed. Each time sample a new one at
% random
for i = 1:nets
        
    %check time
    sprintf(' Finished Net %d\n', i)
    cputime-ttt
    ttt=cputime;
    
    % randomize the bayesian net
    bnet = create_random_bnet(nodes, ep);  % How do we randomize ????? 
    
    % Find the edge distances of all the dags
    if(do_all_dags)        
        for j=1:length(exh_dags)
            edge_all_dags_exh_score(i, j) =  dg_dag_dist(exh_dags{j}, bnet.dag, undirected);
        end


    % New : Find which DAGs are contained in the DAG of the BNET generating
    % the data


    end
        
    number_of_dags_now_is = length(exh_dags)
    
    % Save the index of the correct dag
    for j=1:length(exh_dags)
        if(  sum(sum(bnet.dag == exh_dags{j})) == nodes*nodes  )
            correct_index = j;
            break;
        end
    end
    
    % Here we enumerated all dags and didn't find exh_dags
    if(j == length(exh_dags))
       if( sum(sum(bnet.dag == exh_dags{j})) ~= nodes*nodes )
           correct_index = -1;
           sprintf('Error, couldnt find the DAG');
       end
   end
   
    sprintf('The correct DAG is index %d\n', correct_index); % How come that sometimes there is no correct index ??? 
  

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
    
    
    sample_ttt=cputime;
    % sample data from it
    for j = 1:max_nsamples
        sam(:,j) = sample_bnet(bnet);
    end
    zuk_sampling_time = cputime-sample_ttt
        
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
   
       % Do exhaustive search for learning ... 
       if(nodes < 6)   
           if(do_exh)    
               

                score_dags_ttt=cputime;
                %bayesian_exh_scores = score_dags(cur_data, 2*ones(1,nodes), exh_dags);
                bayesian_exh_scores = zuk_score_dags(cur_data, 2*ones(1,nodes), exh_dags);  % A faster version
                
                 zuk_score_dags_time = cputime-score_dags_ttt
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
                
%                 res
%                 cur_samples
%                 should_be_int_fraction = cur_samples/res
%                 i
%                 correct_index
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
                
           end   % if do exhaust
       end    % if nodes < 6
       
   
        if(do_mc)
            % Now try to learn ....
            
            % sample maximal number of samples. We will take each time only
            % a part of them ...
            [sampled_graphs, accept_ratio, num_edges] = learn_struct_mcmc(cur_data, 2*ones(1,nodes), 'nsamples', max_ngraph_samples, 'burnin', 10);

   
            % sprintf('In original graph %d edges, in sampled graphs %f edges on average!! \n', sum(sum(bnet.dag)), mean(num_edges)); 
  
            
            for(cur_ngraph_samples = (res:res:max_ngraph_samples))
                
   
                sampled_scores = zuk_score_dags(cur_data, 2*ones(1, nodes), sampled_graphs(1:cur_ngraph_samples) );      
                max_index = find(sampled_scores == max(sampled_scores));
   
                qu = sampled_graphs(1:cur_ngraph_samples);
 
   
                % now calculate an 'average' DAG
                ave_dag = zeros(nodes);
                sg=zeros(nodes,nodes,cur_ngraph_samples);
                for k=1:cur_ngraph_samples
                %%%      sg(:,:,k) = cell2num(sampled_graphs(1)); sprintf('SG %d\n', sg(:,:,k) ) 
                    sg(:,:,k) = sampled_graphs{k};
                end   
   
   
                sgu = reshape(sg, nodes*nodes, cur_ngraph_samples);
   
                % Now do the unique 
                sgu = unique(sgu', 'rows')';
   
   
                sgu = reshape(sgu, nodes, nodes, size(sgu,2));
   
                for k=1:size(sgu,3)
                    dag_u{k} = sgu(:,:,k);
                end

                mcmc_post = mcmc_sample_to_hist(sampled_graphs(1:cur_ngraph_samples), dag_u);   

                
                % Do a new ave !!!!!!
                new_ave_dag = zeros(nodes);
                for k=1:length(mcmc_post)
                    new_ave_dag = new_ave_dag + (dag_u{k} * mcmc_post(k));
                end
   
       
                new_ave_dag = new_ave_dag / sum(mcmc_post);
                new_ave_dag(new_ave_dag > 0.5) = 1;
                new_ave_dag = floor(new_ave_dag);

   
                ave_dag(sum(sg,3) > 0.5) = 1;
    
    
                max_learn_score(cur_samples/res, cur_ngraph_samples/res) = max_learn_score(cur_samples/res, cur_ngraph_samples/res) + ...
                dg_dag_dist(qu{max_index(1)}, bnet.dag, undirected);
                ave_learn_score(cur_samples/res, cur_ngraph_samples/res) = ave_learn_score(cur_samples/res, cur_ngraph_samples/res) + ...
                dg_dag_dist(ave_dag, bnet.dag, undirected);
                new_ave_score(cur_samples/res, cur_ngraph_samples/res) = new_ave_score(cur_samples/res, cur_ngraph_samples/res) + ...
                dg_dag_dist(new_ave_dag, bnet.dag, undirected);
 
    
                if(draw_fig)
          
                    figure; draw_graph(qu{max_index(1)}); title('best learned');
                    figure; draw_graph(new_ave_dag); title('new ave learned');
    
                end
            end % enumerating cur_ngraph_samples
        end  % do_mc
        
    end % loop over nsamples
    
%    sprintf('Learned Score : %d\n', sum(sum(qu{max_index(1)} == bnet.dag)))
%    sprintf('AveLear Score : %d\n', sum(sum(ave_dag == bnet.dag)))
   
   
end % loop over nets 



corr_bayesian_edge = zeros(nets, max_nsamples/res);
corr_bayesian_KL = zeros(nets, max_nsamples/res);
% Calculate correlations
if(do_all_dags)
    for i = 1:nets
        for j = 1:max_nsamples/res
        %    sprintf('The Sizes (edge : bayes) \n')
        %    size(edge_all_dags_exh_score)
        %    size(bayesian_all_dags_exh_score)
            
            temp = corrcoef(edge_all_dags_exh_score(i,:), bayesian_all_dags_exh_score(i, j, :));
            corr_bayesian_edge(i, j) = temp(2);
        
            if(do_KL)
         %       sprintf('The Sizes (KL : bayes)\n')
         %       size(KL_all_dags_exh_score)
         %       size(bayesian_all_dags_exh_score)
                temp = corrcoef(KL_all_dags_exh_score(i, j, :), bayesian_all_dags_exh_score(i, j, :));
                corr_bayesian_KL(i, j) = temp(2);
            end
        end
    end

% % %     edge_all_dags_exh_score  = zeros(nets, length(exh_dags)); % Here we save the samples dimention ...
% % %     bayesian_all_dags_exh_score = zeros(nets, max_nsamples/res, length(exh_dags)); 
% % %     KL_all_dags_exh_score  = zeros(nets, max_nsamples/res, length(exh_dags));
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
max_learn_score = max_learn_score/nets; 
ave_learn_score = ave_learn_score/nets;
new_ave_score = new_ave_score/nets;
%edge_max_exh_score = edge_max_exh_score/nets;

KL_max_learn_score = KL_max_learn_score/nets; 
KL_ave_learn_score = KL_ave_learn_score/nets;
KL_new_ave_score = KL_new_ave_score/nets;
%KL_max_exh_score = KL_max_exh_score/nets;




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


% prepare the punish tab. Punishment table for mis-learned edges.
%pun_tab = zeros(3, 3, 3, 3);


          % 00 02 20 11 
pun_tab = [ 0, 3, 3, 3; ... % 00 
            3, 0, 2, 1; ... % 02
            3, 2, 0, 1; ... % 20 
            3, 1, 1, 0];   % 11
        

% Now calculate the score . Punishments for mis-learned edges are : 
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

