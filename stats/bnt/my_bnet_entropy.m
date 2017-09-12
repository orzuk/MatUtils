function bn_ent = my_bnet_entropy(bnet)
% my_BNET_ENTROPY computes the entropy of a distribution given by a bayesian network, 
%
%     ---
%     \
%     /    P(x) * log(P(x))
%     ---
%      x
%
%
%  Clearly, all nodes are discrete. For now we support only boolean nodes.
%
%

epsilon = 0.000000000001;


% First calculate the entropy of bnet1
% .....
nodes = length(bnet.dag);

bn_ent = 0;



% Choose inference engine. Note : this enumerative engine is very slow, and
% should be replaced by a different one. I tried to use jtree or var_elim,
% but for jtree there are sets of nodes which cannot be margenalized, and
% for var_elim, which seemed to me the best engine to use, i got an error
% (something about no T2 variable) - could be a bug!
engine = enumerative_inf_engine(bnet);
evidence = cell(1, nodes);
engine = enter_evidence(engine, evidence);



% Open all the nodes ..
m = marginal_nodes(engine, [1:nodes]);
probs = (reshape(m.T,  prod(bnet.node_sizes), 1) + epsilon) / (1.0 + epsilon);

debug_entropy = entropy(probs);
bn_ent = debug_entropy;



% % for i=1:nodes
% %         
% %     cur_ent = 0;
% %     
% %     % Here we compute the conditional entropy of x_i, given all the
% %     % possibilities for its parents
% %     ps = parents(bnet.dag, j);
% %     
% %     % Calculate the probability distribution of the parents using the engine 
% %     m = marginal_nodes(engine, ps);
% %     
% %     m.T;
% %     
% %     
% %     % i is a root
% %     if(isempty(ps))
% %         for j=1:ns(i)
% %             cur_ent = cur_ent + bnet.CPD{i}(j) * log( bnet.CPD{i}(j) );
% %         end
% %     end
% %               
% %     
% % 
% % end






