% Find the DAG (model) which is closest in relative entropy to
% a given distribution (point). We allow the maximum in-degree
% to be max_k. We currently check only DAGs with the  'natural' orientation.
% Compute the conditional mutual information of X and Y given Z,
% Where their joint distribution is P
function [MI] = conditional_mutual_information(P, X,Y,Z)
n = log2(size(P,2)); % Get the dimenstion

H_Z = entropy(collapse_prob(P,Z)');
H_XZ = entropy(collapse_prob(P,union(X,Z))');
H_YZ = entropy(collapse_prob(P,union(Y,Z))');
H_XYZ = entropy(collapse_prob(P,union(union(X,Y),Z))');

MI = H_XZ + H_YZ - H_Z - H_XYZ;





