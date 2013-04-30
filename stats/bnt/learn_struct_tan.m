function dag = learn_struct_tan(data, class_node, root, node_sizes, scoring_fn)
% LEARN_STRUCT_TAN Learn the structure of the tree augmented naive bayesian network 
% (with discrete nodes)
% dag = learn_struct_tan(app, class, root,node_sizes,scoring_fn)
%
% Input :
% 	data(i,m) is the value of node i in case m
% 	class_node is the class node
% 	root is the root node of the tree part of the dag (must be different from the class node)
%   	node_sizes = 1 if gaussian node,
%   	scoring_fn = 'bic' (default value) or 'mutual_info'
%
% Output :
%	dag = adjacency matrix of the dag
%
% V1.1 : 21 may 2003, (O. Francois, Ph. Leray)


if nargin <4
    error('Requires at least 4 arguments.')
end

if nargin == 4
    scoring_fn='bic';
end;

if class_node==root
  error(' The root node can''t be the class node.');
end

if root>class_node
  root=root-1;
end

N=size(data,1);
node_type=cell(N-1,1);
for i=1:N-1
  node_type{i}='tabular';
end

dag=zeros(N);
notclass=setdiff(1:N,class_node);
T = learn_struct_mwst(data(notclass,:), ones(1,N-1), node_sizes(notclass), node_type, scoring_fn, root);
dag(class_node,notclass)=1;
dag(notclass,notclass)=T;