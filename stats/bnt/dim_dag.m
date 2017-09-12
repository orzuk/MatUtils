function dim = dim_dag( G)
% Get the dimension of a binary dag
dim=sum(2.^sum(G));

