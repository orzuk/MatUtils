% Compute for many (discrete noded) bnets, their joint probability
% distribution. Currently works only for binary nodes
function bnets_dists = my_bnet_to_probs(bnets)
% Choose inference engine. Note : this enumerative engine is very slow, and
% should be replaced by a different one. I tried to use jtree or var_elim,
% but for jtree there are sets of nodes which cannot be margenalized, and
% for var_elim, which seemed to me the best engine to use, i got an error
% (something about no T2 variable) - could be a bug!
nodes = length(bnets{1}.dag);
evidence = cell(1, nodes);
epsilon = 0.000000001;
bnets_dists = zeros(length(bnets), 2^nodes);

% Open all the nodes ..
for j=1:length(bnets)
%%    engine = enumerative_inf_engine(bnets{j});
    engine = var_elim_inf_engine(bnets{j});
    engine = enter_evidence(engine, evidence);
    m = marginal_nodes(engine, [1:nodes]);
    bnets_dists(j,:) = reshape(m.T,  prod(bnets{j}.node_sizes), 1); % apply Dirichlet correction
end

% Now normalize the joint probability
bnets_dists = (bnets_dists + epsilon) ./ (1.0 + epsilon);
bnets_dists = bnets_dists ./ repmat(sum(bnets_dists,2), 1, 2^nodes);