function bnet = create_random_bnet(nodes, edge_prob)
% create_random_bnet creates a 'random' bayesian net.
% The number of input nodes is nodes.
% For each pair of nodes (i,j), there will be an edge with probability
% edge_prob. This is : edge_prob/2 for i->j and edge_prob/2 for j->i
% Method : 
% 1. Create a random DAG
% 2. Create random CPD's, with the standard derichlet distribution on the parameters
%

% create random DAG
dag = create_random_DAG(nodes, edge_prob);


% Now make a bayesian networks from the dag. 
% For now, we assume all nodes are discrete and boolean
ns = 2*ones(1,nodes); 

bnet = mk_bnet(dag, ns);


% Now account for the CPD's
 derich_p = 1.0;
    
% randomize the CPD's
for j=1:nodes
    ps = parents(bnet.dag, j);
        
    psz = power(2, length(ps)); %prod(ns(ps));
    CPT = dirichlet_sample(derich_p*ones(1,2), psz)
    bnet.CPD{j} = tabular_CPD(bnet, j, 'CPT', CPT);
    
end
