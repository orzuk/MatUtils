function bnets = create_bounded_degree_bnet(nodes, max_deg, num_bnets)
% create_random_bnet creates a 'random' bayesian net.
% The number of input nodes is nodes.
% For each pair of nodes (i,j), there will be an edge with probability
% edge_prob. This is : edge_prob/2 for i->j and edge_prob/2 for j->i
% Method :
% 1. Create a random DAG
% 2. Create random CPD's, with the standard derichlet distribution on the parameters
%

% create random DAG
dags = create_bounded_degree_DAG(nodes, max_deg, num_bnets);

% %  Add a random edge
% for i=1:num_bnets
%     Mi = min(find(dags(1,2:end,i) == 0));
%     if(~isempty(Mi))
%         dags(1,Mi+1,i) = 1;
%     end
% end

% Now make a bayesian networks from the dag.
% For now, we assume all nodes are discrete and boolean
ns = 2*ones(1,nodes);

% Now account for the CPD's
derich_p = 1.0;

% randomize the CPD's
bnets = [];
for i=1:num_bnets
    bnets{i} = mk_bnet(dags(:,:,i), ns);
    for j=1:nodes
        ps = parents(bnets{i}.dag, j);
        psz = power(2, length(ps));
        CPT = dirichlet_sample(derich_p*ones(1,2), psz);
        bnets{i}.CPD{j} = tabular_CPD(bnets{i}, j, 'CPT', CPT);
    end
    
    bnets{i}.my_CPT = zeros(nodes,2);
    bnets{i}.my_CPT(bnets{i}.order(1),:) = CPD_to_CPT(bnets{i}.CPD{bnets{i}.order(1)})';
    for j=2:nodes
        tmp = CPD_to_CPT(bnets{i}.CPD{bnets{i}.order(j)})';
        bnets{i}.my_CPT(bnets{i}.order(j),:)  =  tmp(1,:);
    end

end



% % % % % Alternative: do it ourselves: (works only for binary nodes)
% % % % bnets = [];
% % % % bnets.dags = dags;
% % % % bnets.dims = reshape(sum(2.^sum(dags)), 1, num_bnets);
% % % % bnets.CPD = rand(bnets.dims(1),num_bnets); 
% % % % 
% % % % 
% % % % % Get already the probs, 'on the way'
% % % % P = ones(2^nodes,num_bnets);
% % % % % Try with a loop first
% % % % 
% % % % parents = find(dags
% % % % 
% % % % BitsVec = zeros(2^nodes, nodes);
% % % % for i=1:nodes
% % % %     BitsVec(:,i) = bitget(0:2^nodes-1,i);
% % % % end    
% % % % 
% % % % for i=1:num_bnets
% % % %     for j=1:nodes
% % % %        
% % % %         
% % % %     end
% % % %     
% % % % end
% % % % 
% % % % 
% % % % 
% % % % 
