function children = get_children( X, G)
% Get the children of X in G
[I J] = find(G(X,:)); 
children = setdiff(J, X); 
% parents = find(G(X,:)); 
return;
