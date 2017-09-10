% Calculate the distance matrix between two sets of vectors
function D = two_vecs_to_distance(V1, V2)

D = repmat(sum(V1.^2,2), 1, size(V2, 1)) + repmat(sum(V2.^2,2), 1, size(V1, 1))';
D = sqrt( max(0, D - 2 .* V1 * V2') ); % rounding errors might give negatives, for singles 
%% D = sqrt( D + D' - 2 .* V * V' ); % rounding errors might give negatives


