%%%function [max_learn_score, ave_learn_score, new_ave_score, max_exh_score, KL_max_learn_score, KL_ave_learn_score, KL_new_ave_score, KL_max_exh_score, correct_rank] = ...
%%% New function ! 28.2.06 (Zuk)
function [edge_max_exh_score, KL_max_exh_score, edge_average_exh_score, correct_rank, edge_all_dags_exh_score, ...
    bayesian_all_dags_exh_score, KL_all_dags_exh_score, corr_bayesian_edge, corr_bayesian_KL, exh_dags, ...
    dags_order_mat, containing_error, not_containing_error] = ...
    check_learning_quality(bnet, nodes, iters, max_nsamples, max_ngraph_samples, res, do_exh, do_mc, do_all_dags, do_KL, draw_fig)
% create random samples from a GIVEN bayesian nets, and then try to learn the net.
% We want to see how well we learn certain things
% The number of input nodes is nodes.
% max_nsamples is the maximal number of samples we draw
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
TOL= 0.00000000000001;

nets = 1; % Artificial, to be removed

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
    num_dags = length(exh_dags)
end



% New : this vector says if the true probability distribution is contained
% in the DAG
iss_contained_flags = zeros(1, length(exh_dags));


ttt=cputime;

% Generate samples from our bayesian net
sprintf(' Starting Net')
cputime-ttt
ttt=cputime;

% Find the edge distances of all the dags
if(do_all_dags)
    for j=1:length(exh_dags)

        %         cur_dag_j = exh_dags{j}
        %         bnet_dag = bnet.dag

        edge_all_dags_exh_score(1, j) =  dg_dag_dist(exh_dags{j}, bnet.dag, undirected);
    end






    % At the beginning find all the equivalence classes
    equiv_class_indexes = zeros(1,num_dags); % Save enough space
    equiv_class_vec = zeros(1,num_dags)-1;
    cur_eq=0;
    for j=1:num_dags

        if(equiv_class_vec(j) == (-1))
            equiv_class_vec(j)=cur_eq; cur_eq=cur_eq+1; equiv_class_indexes(cur_eq) = j; % Save also the first of each class
        end
        for k=j+1:num_dags
            if(are_dags_equiv( exh_dags{j}, exh_dags{k}))
                if((equiv_class_vec(k) == (-1)) )
                    equiv_class_vec(k) = equiv_class_vec(j);
                end
            end
        end
    end

    num_equiv_classes = cur_eq
    equiv_class_indexes = equiv_class_indexes(1:num_equiv_classes); % Save only the number of equivalence classes.

    % Shrink everything to equivalence classes
    exh_dags = exh_dags(equiv_class_indexes); save_exh_dags = exh_dags;

    % New : Find which DAGs are contained in the DAG of the BNET generating
    % the data
    dags_order_mat = zeros(num_equiv_classes);

    dag_indep_mat = zeros(num_equiv_classes, 4^nodes);
    for j=1:num_equiv_classes
        dag_indep_mat(j,:) = get_dag_indep(exh_dags{j}); % Find all independences from a given DAGs.

        if(mod(j,50) == 0)
            doing_dag_j = j
        end
    end


    for j=1:num_equiv_classes
        for k=1:num_equiv_classes

            % Check if j-th graph contains the k-th graph
            if( isempty(find(dag_indep_mat(j,:) < dag_indep_mat(k,:))) )
                dags_order_mat(j,k) = 1;
            end
        end
    end



    % Find which DAGs is the P-map for the bnet generating the data.
    my_dag_index=-1
    for j=1:num_equiv_classes
        if(are_dags_equiv(bnet.dag, exh_dags{j}))
            my_dag_index = j
        end
    end
    correct_index=my_dag_index;

    % Now count how many dags are contained in this DAG
    num_DAGS_contained = sum(dags_order_mat(my_dag_index,:))
    num_DAGS_contained_reverse = sum(dags_order_mat(:,my_dag_index))






    % Seperate the dags (actually equivalence classes) according to the
    % relations with the correct dag
    larger_dags = find(dags_order_mat(my_dag_index,:))'
    smaller_dags = find(dags_order_mat(:,my_dag_index))
    equiv_dags = intersect(larger_dags, smaller_dags)
    %%%% equiv_dags = find(equiv_class_vec == equiv_class_vec(my_dag_index))
    strictly_larger_dags = setdiff(larger_dags, equiv_dags);
    strictly_smaller_dags = setdiff(smaller_dags, equiv_dags);


    containing_dags = strictly_larger_dags;
    not_containing_dags = setdiff(1:num_equiv_classes, containing_dags);
    not_containing_dags = setdiff(not_containing_dags, equiv_dags);


    % Check symmetry
    is_zero_min = min(min(dags_order_mat - dags_order_mat'))
    is_zero_max = max(max(dags_order_mat - dags_order_mat'))

    size(larger_dags)
    size(smaller_dags)
    size(equiv_dags)



    % Do a quick plot
    figure; subplot(2,2,1); hold on; graph_draw(exh_dags{strictly_larger_dags(1)}); title('Larger');
    subplot(2,2,2); hold on; graph_draw(exh_dags{strictly_smaller_dags(1)}); title('Smaller');
    subplot(2,2,3); hold on; graph_draw(exh_dags{equiv_dags(1)}); title('Equivalent');
    subplot(2,2,4); hold on; graph_draw(bnet.dag); title('Original');



    % Here take only one smaller DAG and one larger DAG, to save time
    take_only_one_flag=1;

    if(take_only_one_flag)

        % Take minimal dimension of containing dags
        best_j=-1; best_dim=2^nodes+1;
        strictly_larger_dags
        for j=strictly_larger_dags'
            %             j
            %             exh_G = exh_dags{j}
            if( dim_dag(exh_dags{j}) < best_dim)
                best_j=j; best_dim=dim_dag(exh_dags{j});
            end
        end
        larger_dag_index = best_j;


        % Take maximal dimension of smaller dags
        best_j=-1; best_dim=-1;
        for j=strictly_smaller_dags'
            if( dim_dag(exh_dags{j}) > best_dim)
                best_j=j; best_dim=dim_dag(exh_dags{j});
            end
        end
        smaller_dag_index = best_j;

        % Shrink to get only 3 dags: The correct one, larger and smaller
        length(exh_dags)
        my_dag_index
        larger_dag_index
        smaller_dag_index
        exh_dags = exh_dags([my_dag_index, larger_dag_index, smaller_dag_index]);
        correct_index = 1; containing_dags = 2; not_containing_dags = 3; num_equiv_classes=3;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform special pre-computations for the 3-graphs in question:
        Q  = {}; % Holds the distributions of the 3 bnets. Good only for small binary nets.

        % compute the projection of P on  the smaller DAG which is
        Q{1} = zuk_bnet_to_probs(bnet); % This is the correct distribution
        Q{2} = Q{1};  % Containing graph also contains the correct distribution
        Q{3} = project_on_bnet(Q{1}, exh_dags{3});


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    % Here loop over all Dags and ... Default for do_all_dags is one !
    edge_all_dags_exh_score  = zeros(nets, num_equiv_classes); % Here we save the samples dimention ...
    bayesian_all_dags_exh_score = zeros(nets, max_nsamples/res, num_equiv_classes);
    KL_all_dags_exh_score  = zeros(nets, max_nsamples/res, num_equiv_classes);


end % if do_all_dags


number_of_dags_now_is = length(exh_dags)



sam = cell(nodes, max_nsamples);


sample_ttt=cputime;




if(draw_fig)
    figure; draw_graph(bnet.dag); title('correct');
    figure; draw_graph(dag_to_cpdag(bnet.dag)); title('correct cpdag');
end




nsamples_vec = res:res:max_nsamples;
containing_error = zeros(1,length(nsamples_vec)); % Overfitting: A larger graph containing the true dag was learned
not_containing_error = zeros(1,length(nsamples_vec)); % Underfitting: A graph not containing the true dag was learned
overall_error = zeros(1,length(nsamples_vec)); % Either one of these errors



% Prepare the sampling once for each number of samples
chunk_size = res; % Size of chunk to use.

binom_cumsum_vecs = {}; % Cell array containing what we need for sampling
binom_probs = {};


%for(i = 1:length(nsamples_vec))
%    [bnet_cum_dist, binom_probs, binom_cumsum_vecs{i}] = zuk_sample_bnet_prepare(zuk_bnet_to_probs(bnet), nsamples_vec(i), chunk_size);
%end

% Do preperation for importance sampling
importance_num_dists =10; % Take many importance points, hopefully one of which will give a good sampling
projected_Qs = {}; P=zuk_bnet_to_probs(bnet);
importance_Qs = {};


%%%% Reminder: correct_index = 1; containing_dags = 2; not_containing_dags = 3; num_equiv_classes=3;
if(take_only_one_flag) % Here we need to prepare for doing importance sampling
    for g=1:length(exh_dags)
        projected_Qs{g} = project_on_bnet( P, exh_dags{g});

        ggg=g
        QQQ=projected_Qs{g};

        % Do correction for the containing graph !!!!  Take a random distribution !!!
        if(g==2)
            R = rand(length(P),1); R = R ./ sum(R);
            projected_Qs{g} = project_on_bnet( R, exh_dags{g});
        end

        for q=1:importance_num_dists

            % Compute the 'middle' distribution - take weighted average
            importance_Qs{g,q} = ( P*q + projected_Qs{g}*(importance_num_dists-q) ) / importance_num_dists;

            for(i = 1:length(nsamples_vec))
                [bnet_cum_dist, binom_probs{g,q,i}, binom_cumsum_vecs{g,q,i}] = ...
                    zuk_sample_bnet_prepare(importance_Qs{g,q}, nsamples_vec(i), chunk_size);
            end


        end
    end
end


% Run the learning many times to get average probability
score_dags_ttt=cputime;


% This is already in size iters*2^n

data_counts = {};

sample_ttt=cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Here sample directly each point
% Skip the original DAG to save time!

if(take_only_one_flag)
    run_over_dags_vec = 2:length(exh_dags);

    for g=run_over_dags_vec % For each graph we do a different importance sampling suited for him!

        for q=1:importance_num_dists

            data_counts{g,q,1} =  zuk_sample_bnet_generate_many_samples(bnet_cum_dist, res, iters, binom_probs{g,q,1},binom_cumsum_vecs{g,q,1});
            for i=2:length(nsamples_vec)
                data_counts{g,q,i} =  data_counts{g,q,i-1} + ...
                    zuk_sample_bnet_generate_many_samples(bnet_cum_dist, res, iters, binom_probs{g,q,i},binom_cumsum_vecs{g,q,i});
            end
        end
    end

else
    run_over_dags_vec = 1:length(exh_dags);

    for(i = 1:length(nsamples_vec))
        [bnet_cum_dist, binom_probs{1,1,i}, binom_cumsum_vecs{1,1,i}] = ...
            zuk_sample_bnet_prepare(P, nsamples_vec(i), chunk_size);
    end

    % Do only one count , with the original Q !!!!!!!!
    data_counts{1,1,1} =  zuk_sample_bnet_generate_many_samples(bnet_cum_dist, res, iters, binom_probs{1,1,1},binom_cumsum_vecs{1,1,1});
    for i=2:length(nsamples_vec)
        data_counts{1,1,i} =  data_counts{1,1,i-1} + ...
            zuk_sample_bnet_generate_many_samples(bnet_cum_dist, res, iters, binom_probs{1,1,i},binom_cumsum_vecs{1,1,i});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Here sample using mcmc
%
% num_mix_steps = 500;
% for g=1:length(exh_dags) % For each graph we do a different importance sampling suited for him!
%     for i=1:length(nsamples_vec)
%         % Generate one init point
%         init_count = zuk_sample_bnet_generate_many_samples(bnet_cum_dist, nsamples_vec(i), 1, binom_probs,binom_cumsum_vecs{g,i});
%         data_counts{g,i} = zuk_mcmc_sample_bnet_counts(init_count, num_mix_steps, bnet_cum_dist, iters, binom_cumsum_vecs);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



all_zuk_sampling_time = cputime-sample_ttt

% size_data_counts = size(data_counts{1})
% size_data_counts2 = size(data_counts{length(nsamples_vec)})
%
% zuk_score_dags_time = cputime-score_dags_ttt;
score_dags_ttt=cputime;

containing_error = zeros(1, length(nsamples_vec));
not_containing_error = zeros(1, length(nsamples_vec));

% Go and increase the samples, each time by res
for(j = 1:length(nsamples_vec))

    cur_samples = nsamples_vec(j)
    cur_index = cur_samples/res;

    % Work only on part of the data

    % Do exhaustive search for learning ...
    if(nodes < 6)
        if(do_exh)

            %%%% Reminder: correct_index = 1; containing_dags = 2; not_containing_dags = 3; num_equiv_classes=3;

            if(take_only_one_flag) % Work only with one containing and one not-containing graphs !!!!

                for q=1:importance_num_dists  % Loop over different importance distributions

                    bayesian_exh_scores = zuk_BIC_score_dags(data_counts{2,q,j}, nsamples_vec(j), exh_dags);% A faster MDL version
                    importance_correction =zeros(iters,1);

                    for i=1:length(P)
                        importance_correction = importance_correction + ( log(P(i)) - log(importance_Qs{2,q}(i)) ) .* data_counts{2,q,j}(:,i);
                    end
                    importance_correction = exp(importance_correction);
                    size(importance_correction)
                    size( bayesian_exh_scores(:,2))
                    containing_error(j) = containing_error(j) + sum( (bayesian_exh_scores(:,1) < bayesian_exh_scores(:,2)) .* importance_correction );


                    containing_not_important_prob = sum( (bayesian_exh_scores(:,1) < bayesian_exh_scores(:,2)) ) / iters;


                    bayesian_exh_scores = zuk_BIC_score_dags(data_counts{3,q,j}, nsamples_vec(j), exh_dags);  % A faster MDL version
                    importance_correction = zeros(iters,1);
                    for i=1:length(P)
                        importance_correction = importance_correction + ( log(P(i)) - log(importance_Qs{3,q}(i)) ) .* data_counts{3,q,j}(:,i);
                    end
                    importance_correction = exp(importance_correction);
                    not_containing_error(j) = not_containing_error(j) + sum( (bayesian_exh_scores(:,1) < bayesian_exh_scores(:,3)) .* importance_correction);

                    not_containing_not_important_prob = sum( (bayesian_exh_scores(:,1) < bayesian_exh_scores(:,3)) ) / iters;


                end

            else % Here work with ALL equivalence classes !!!!

                % A faster MDL version - here all DAGS! and we use the
                % counts from the true Q !!!

                two_ttt = cputime;
                
                bayesian_exh_scores = zuk_BIC_score_dags(data_counts{1,1,j}, nsamples_vec(j), exh_dags);
                
%                 ALL_scores_time = cputime - two_ttt
%                 
%                 bayesian_exh_scores = zuk_BIC_score_dags(data_counts{1,1,j}, nsamples_vec(j), exh_dags(20))
%                 
%                 OLD_score_time = cputime - two_ttt - ALL_scores_time
%                 
%                 size(bayesian_exh_scores)
%                 size(bayesian_exh_scores2)
%                 
%                 
%                bayesian_exh_scores2-bayesian_exh_scores
%                
%                MMMM = max(max((bayesian_exh_scores2-bayesian_exh_scores) ./ bayesian_exh_scores))
%                mmmm = min(min((bayesian_exh_scores2-bayesian_exh_scores) ./ bayesian_exh_scores))
%                 
%                 mean_correct = mean(bayesian_exh_scores)
%                 mean_farsh = mean(bayesian_exh_scores2)
%                 maxmax = max(max(bayesian_exh_scores-bayesian_exh_scores2))
%                 minmin = min(min(bayesian_exh_scores-bayesian_exh_scores2))
                %                 size(bayesian_exh_scores(:,my_dag_index))
                %                 size(bayesian_exh_scores(containing_dags,:))
                %                 size( max(bayesian_exh_scores(containing_dags,:)) )
                %                 size( max(bayesian_exh_scores(containing_dags,:), [],  2) )
                %

                max_exh_scores = max(bayesian_exh_scores);
                
               overall_error(j) = sum(bayesian_exh_scores(my_dag_index,:) < max(bayesian_exh_scores)  );
%                 containing_error(j) = sum(max(bayesian_exh_scores(containing_dags,:)) < max_exh_scores );
%                 not_containing_error(j) = sum(max(bayesian_exh_scores(not_containing_dags,:)) < max_exh_scores );
                containing_error(j) = sum(bayesian_exh_scores(my_dag_index,:) < max(bayesian_exh_scores(containing_dags,:)) );
                not_containing_error(j) = sum(bayesian_exh_scores(my_dag_index,:) < max(bayesian_exh_scores(not_containing_dags,:)) );
                
                
            end


        end   % if do exhaust
    end    % if nodes < 6

end % loop over nsamples


all_zuk_score_dags_time = cputime-score_dags_ttt

% Normalize
if(take_only_one_flag) %
    containing_error = containing_error./(iters*importance_num_dists);
    not_containing_error = not_containing_error./(iters*importance_num_dists);
else
    containing_error = max(containing_error./iters , TOL);
    not_containing_error = max(not_containing_error./iters , TOL);
    overall_error = max(overall_error./iters, TOL);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% New figure: Combining error
% figure; hold on;
% plot(nsamples_vec, 1-containing_error-not_containing_error, 'r');
% plot(nsamples_vec, 1-not_containing_error);
% xlabel('num samples'); ylabel('Prob.'); legend('true', 'true+overfit');

figure; hold on;
area(nsamples_vec, [1-overall_error; containing_error; not_containing_error]');
xlabel('num samples'); ylabel('Prob.');  legend('G^*', 'non-minimal I-maps', 'Not I-maps');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Before doing figure, do fit to get the constants missing in Haughton's paper

% SemiLog fit
containing_error
not_containing_error
size(containing_error)
size(nsamples_vec)
[fit_not_containing qual_not_containing] = fit(nsamples_vec', log(containing_error'), 'poly1'); % alpha in the exponent

% LogLog fit
[fit_containing qual_containing] = fit(log(nsamples_vec'), log(containing_error'), 'poly1');


% Here try to calculate what we think will be the constants, from the paper:

% Calc the constant in the algebraic error of overfitting
woodroofe_const = (dim_dag(exh_dags{1}) - dim_dag(exh_dags{2}))/2
fit_woodroofe_const = -fit_containing.p1


% Calc the relative entropy of the true distribution from the underfit ..

if(take_only_one_flag) %
    haughton_underfit_const_opposite = relative_entropy(P, projected_Qs{3}) * log(2)
    haughton_underfit_const = relative_entropy(projected_Qs{3}, P) * log(2) % Note: This is only an UPPERBOUND !Perhaps there is a closer dist. to P
end

fit_haughton_underfit_const = -fit_not_containing.p1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% New figure: Divide error into two cases: containing and not
%%%%% containing dags
figure; hold on; subplot(2,2,1); hold on; plot(nsamples_vec, containing_error, '*');
plot(nsamples_vec, not_containing_error, 'r+');xlabel('num samples'); ylabel('Pr(error)');
title('Two error probabilities as a function of sample size');
legend('Larger', 'Smaller');




% plot the same with my loglog
subplot(2,2,2);  hold on; plot(nsamples_vec, log(containing_error), '*');
plot(nsamples_vec, log(not_containing_error), 'r+');xlabel('num samples'); ylabel('Log Pr(error)');
title('SemiLogs of Two error probabilities as a function of sample size');
legend('Larger', 'Smaller');




% plot the same with loglog
subplot(2,2,3); loglog(nsamples_vec, containing_error, '*'); xlabel('num samples'); ylabel('Pr(error)');
title('Log Larger error probabilities as a function of sample size');
legend('Larger');


subplot(2,2,4); loglog(nsamples_vec, not_containing_error, 'r+'); xlabel('num samples'); ylabel('Pr(error)');
title('Log Smaller error probabilities as a function of sample size');
legend('Smaller');

%%%%% End of new figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


corr_bayesian_edge = zeros(nets, max_nsamples/res);
corr_bayesian_KL = zeros(nets, max_nsamples/res);
% Calculate correlations
% % % % % if(do_all_dags)
% % % % %     for k = 1:nets
% % % % %         for j = 1:max_nsamples/res
% % % % %
% % % % %             temp = corrcoef(edge_all_dags_exh_score(k,:), bayesian_all_dags_exh_score(k, j, :));
% % % % %             corr_bayesian_edge(k, j) = temp(2);
% % % % %
% % % % %             if(do_KL)
% % % % %                 temp = corrcoef(KL_all_dags_exh_score(k, j, :), bayesian_all_dags_exh_score(k, j, :));
% % % % %                 corr_bayesian_KL(k, j) = temp(2);
% % % % %             end
% % % % %         end
% % % % %     end
% % % % %
% % % % %     % % %     edge_all_dags_exh_score  = zeros(nets, length(exh_dags)); % Here we save the samples dimention ...
% % % % %     % % %     bayesian_all_dags_exh_score = zeros(nets, max_nsamples/res, length(exh_dags));
% % % % %     % % %     KL_all_dags_exh_score  = zeros(nets, max_nsamples/res, length(exh_dags));
% % % % % end



if(nodes < 6)
    if(do_exh)
        sprintf('Exhaustive Search Score : %d\n', edge_max_exh_score);
    end
end


% Here we return things
max_learn_score = max_learn_score/nets;
ave_learn_score = ave_learn_score/nets;
new_ave_score = new_ave_score/nets;
%edge_max_exh_score = edge_max_exh_score/nets;

KL_max_learn_score = KL_max_learn_score/nets;
KL_ave_learn_score = KL_ave_learn_score/nets;
KL_new_ave_score = KL_new_ave_score/nets;
%KL_max_exh_score = KL_max_exh_score/nets;


% Return all equiv. classes
exh_dags = save_exh_dags;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare two DGs to see how many edges differ. The first is the directed
% graph we've got. The second is a DAG - usually the original DAG that we
% sampled from. If directed is on, we account also for the direction of
% edges learned. Otherwise we account only for edges/non edges
function dd = dg_dag_dist(dg1, dag2, directed)


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

