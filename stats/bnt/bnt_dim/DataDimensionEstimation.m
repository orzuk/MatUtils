% Compute the intrinsic dimension of a data.
% The data is of size n*m and is assumed to lie on a lower-dimension manifold 
% which the function tries to estimate 
function d = EstimateDimension(D)

metric = 0; % euclidian
knn = size(D, 1); % number of points
[KNN_sim KNN_inds] = KNNBruteForce(D, [], metric, knn); 

