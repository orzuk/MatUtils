function dags = create_bounded_degree_DAG(nodes, max_deg, num_dags)
% create_bounded_degree_DAG creates 'random' directed acyclic graphs.
% The number of input nodes is nodes.
% The maximum (in)-degree of the nodes is max_deg.


dags = zeros(nodes,nodes,num_dags);

% Take all possible sources at the beginning
for i=1:max_deg+1
    dags([1:i-1],i,:) = 1;
end
% randomize sources from now on
for i=max_deg+2:nodes
    for j=1:num_dags
        P = randperm(i-1);
        dags(P(1:max_deg),i,j) = 1;
    end
end
    
% % Permute or not permute?
for j=1:num_dags
    p = randperm(nodes); % create order permutation
    dags(:,:,j) = dags(p,p,j); % permute
end