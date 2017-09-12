function parents = get_parents( X, G)
% Get the parents of X in G
[I J] = find(G(:,X)); 
parents = setdiff(I, X); 
return;
