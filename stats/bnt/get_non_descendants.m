function non_des = get_non_descendants( X, G)
% Get all the non-descendant of X in G
% As a matter of convention, we remove also the parents here!

N = length(G);
non_des = setdiff(1:N,get_descendants( X, G));
non_des = setdiff(setdiff(non_des,X), get_parents(X,G));

return;