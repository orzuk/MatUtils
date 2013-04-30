function rel_ent = my_bnet_relative_entropy(bnet1, bnet2)
% my_BNET_RELATIVE_ENTROPY computes the relative entropy , 
% the Kullbak-Leibner distance between two distributions, defined by:
%
%     ---
%     \
%     /    P(x) * log(P(x)/Q(x)
%     ---
%      x
%
%  Note : This distance is not symetric. Order of arguments is thus
%  important.
%
% G{i} is the i'th dag
%
% Note: the number of DAGs is super-exponential in N, so don't call this with N > 4.


% First calculate the entropy of bnet1
rel_ent = my_bnet_entropy(bnet1);


% Now calculate the other part ...

epsilon = 0.000000000001;


% First calculate the entropy of bnet1
% .....
nodes = length(bnet1.dag);

bn_ent = 0;

% Choose inference engine. Note : this enumerative engine is very slow, and
% should be replaced by a different one. I tried to use jtree or var_elim,
% but for jtree there are sets of nodes which cannot be margenalized, and
% for var_elim, which seemed to me the best engine to use, i got an error
% (something about no T2 variable) - could be a bug!
engine = enumerative_inf_engine(bnet1);
evidence = cell(1, nodes);
engine = enter_evidence(engine, evidence);

% Open all the nodes ..
m = marginal_nodes(engine, [1:nodes]);
Pprobs = (reshape(m.T,  prod(bnet1.node_sizes), 1) + epsilon) / (1.0 + epsilon);

engine = enumerative_inf_engine(bnet2);
evidence = cell(1, nodes);
engine = enter_evidence(engine, evidence);

% Open all the nodes ..
m = marginal_nodes(engine, [1:nodes]);
Qprobs = (reshape(m.T,  prod(bnet2.node_sizes), 1) + epsilon) / (1.0 + epsilon);


debug_entropy = cross_entropy(Pprobs, Qprobs);

rel_ent = debug_entropy; 

    

