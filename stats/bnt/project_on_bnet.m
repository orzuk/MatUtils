% Find the distribution consistent with graph G
% which is 'closest' to P. We simply do the counting
% according to P. Only works for small binary networks
% The output Q is the one maximizing the likelihood of
% a data whose type is P, of all distributions consistent
% with G. It also minimizes the relative-entropy KL(P||Q)
% of all distributions consistent with G.
function Q = project_on_bnet( P, G)

% The current version works only for graphs in the 'natural' order
G_perm  = topological_sort(G);
GG = G(G_perm,G_perm);

n=size(GG,1); % number of nodes
num_dists = size(P,2); % One can project many distributions at once
% Get the indices transofmration
IndVec = 0:2^n-1; IndVec2 = zeros(1,2^n);
for i=1:n
    IndVec2 = IndVec2+ bitget(IndVec, G_perm(i)) .* 2^(i-1);
end

IndVec(IndVec2+1)  = 1:2^n;
PP = P(IndVec,:); % transform the probability P to be according to G's ordering

Q = ones(2^n,num_dists); % Start with ones

% Find the list of relevant vars: i and its parents
cur_vars = {}; cur_vars_inv = {};
cur_vars_bin = zeros(1,n); cur_vars_inv_bin = zeros(1,n); 
for i=1:n % loop over all variables
    cur_vars{i} = [find(GG(:,i))' i];
    cur_vars_inv{i} = setdiff(1:n, cur_vars{i});    
    cur_vars_bin(i) = sum(2.^ (cur_vars{i}-1));
    cur_vars_inv_bin(i) = sum(2.^ (cur_vars_inv{i}-1));
end

% loop over the assignment of all variables - this is the heaviest loop
% (need to replace it by something better !!!)
for x_vec = 0:2^n-1
    for i=1:n % loop over all variables (do we assume here the 'natural' ordering ?
        cur_x_bits = bitget(x_vec, cur_vars{i});
        cur_x_bits_bin = sum(cur_x_bits.* 2.^(0:length(cur_x_bits)-1));
        
        cur_x_bits_inv_x_i = cur_x_bits; cur_x_bits_inv_x_i(end) = 1-cur_x_bits_inv_x_i(end);
        cur_x_bits_inv_x_i_bin = sum(cur_x_bits_inv_x_i.* 2.^(0:length(cur_x_bits_inv_x_i)-1));
        
        cur_indexes = stretch_to_indexes(0:2^length(cur_vars_inv{i})-1, cur_vars_inv_bin(i)) + ...
            stretch_to_indexes(cur_x_bits_bin, cur_vars_bin(i)) + 1;        
        cur_indexes_inv_x_i = stretch_to_indexes(0:2^length(cur_vars_inv{i})-1, cur_vars_inv_bin(i)) + ...
            stretch_to_indexes(cur_x_bits_inv_x_i_bin, cur_vars_bin(i)) + 1;
            
        % loop over all the parents of X_i
        % Do the same loop for many distributions at once
        if(length(cur_indexes) == 1)
            nomin = sum(PP(cur_indexes,:),1);
        else
            nomin = sum(PP(cur_indexes,:));
        end
        if(length(cur_indexes_inv_x_i) == 1)
            denom = nomin + sum(PP(cur_indexes_inv_x_i,:),1);
        else
            denom = nomin + sum(PP(cur_indexes_inv_x_i,:));
        end
%        nomin = sum(PP(cur_indexes,:)) 
%        denom = nomin + sum(PP(cur_indexes_inv_x_i,:))
        Q(x_vec+1,:) = Q(x_vec+1,:) .* nomin ./ denom;
    end
end
% Take the inverse permutation
Q = Q(IndVec2+1,:);


