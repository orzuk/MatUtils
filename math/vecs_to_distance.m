% Calculate the euclidian distance matrix between vectors
%
% Input: 
% V - a matrix mXn of m vectors of size n
% 
% Output: 
% D - a matrix mXm of Euclidian distances between vectors 
% 
function D = vecs_to_distance(V)

D = repmat(sum(V.^2,2), 1, size(V, 1));
D = sqrt( max(0, D + D' - 2 .* V * V') ); % rounding errors might give negatives, for singles 


