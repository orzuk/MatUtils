function is_sub = my_is_submodel(dag1, dag2)
% Check if the equivalence class of dag1 is a submodel of the equivalence
% class of dag2
% Method : For all variables in dag1, check independency in non-ancenstors
% given parents


% Take all Local Markov properties in G2, and check if they are satisfied
% in G2
N = length(dag2);

is_sub = 1; % Assume it is
% Go over all vertices
for i=1:N
    X=i; % The vertex we take
    Y = get_parents(X,dag2); % The parents of X
    Z = get_non_descendants(X,dag2); % The non descendants of X
    if(~dsep(X, Z, Y, dag1))
        is_sub=0;
        return;
    end
end

return;


