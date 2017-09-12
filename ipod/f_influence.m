% Calculate the adoption rate of each vertex given other parameters.
% We calculate lambda given all the constrains : graph, payments etc.
% neighbors_vec is a vector holding the neighbors of each
% vertex
% lambda_0 is the 'basal' rate, i.e. without doing anything
% lambda_1 is the total influence of the neighbors +money+global field
%
function lambda_vec = f_influence(E, X, c_vec, k_global, lambda_0, ...
    lambda_1, money_factor, tv_factor)

N = length(X); % number of vertices
N_on = sum(X); % number of infected vertices

% implementation: take into account only local neighbors,
% money serves as another neighbour and
% global publicity is another neighbour

lambda_vec = zeros(1,N)+lambda_0; % constant lambda_0 rate for everybody

% Compute f in a vectorized way
% X
% E(:,find(X))
% sum(E(:,find(X)))


neighbors = sum(E)--diag(E)' ;
neighbors_on = sum(E(find(X),:),1);


pseudo_neighbours = neighbors + money_factor + tv_factor;
pseudo_on = neighbors_on + tv_factor.*((N_on + k_global) ./ (N + k_global)) + money_factor.*(c_vec./(1+c_vec));

% One possibility : Rate is determined by relative characters of neighbors
%%%%  lambda_vec = lambda_vec + (lambda_1-lambda_0) .* pseudo_on ./ pseudo_neighbours;

% Other possibility : Rate is determined by only absolute contribution as a sum of all neighbors which are on
lambda_vec = lambda_vec + (lambda_1-lambda_0) .* pseudo_on;



% % %   for i=1:N % loop on all vertices. The 'heaviest' part
% % %
% % %     fff = find(E(i,:));
% % %     neighbors = neighbors_vec(i);
% % %     neighbors_on = length(find(X(fff)));
% % %
% % %     pseudo_neighbours = neighbors + money_factor + tv_factor;
% % %     pseudo_on = neighbors_on + tv_factor*((N_on + k_global) / (N + k_global)) + money_factor*(c_vec(i)/(1+c_vec(i)));
% % %
% % % 				% add rate according to pseudo_neighbors
% % %     lambda_vec(i) = lambda_vec(i) + lambda_1 * pseudo_on / pseudo_neighbours;
% % %   end
