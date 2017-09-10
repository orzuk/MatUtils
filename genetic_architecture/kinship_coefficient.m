% Compute kinship coefficient for all pairs of members in a family tree
%
% Input:
% family_tree - tree describing the familial relationships
%
% Output:
% K - matrix of kinship coefficients
% U - matrix of probability of sharing 0,1,2 alleles
%
function [K U] = kinship_coefficient(family_tree)

N = length(family_tree); % number of individuals;

K = eye(N); % set diagonal (self-kinship) to one
U = zeros(N,N,3); % U(i,j,k) = Pr(i and j individuals share k-1 alleles)

ancestors = get_graph_sources(family_tree); % start with ancestors of the family 
U(ancestors, ancestors,1) = 1 - eye(length(ancestors)); % no allele sharing between ancestors
U(:,:,3) = eye(N); % always share 2 alleles with self
U(:,:,2) = 1 - U(:,:,1) - U(:,:,3); % one allele shared

% [member_dists member_discovery_time pred] = bfs(sparse(family_tree), ancestors(1));
%[sorted_member_dists sort_perm] = sort(member_dists); % sort to get ordering (takes care of everything except unreachable nodes)
%first_reachable_member = find(sorted_member_dists > 0, 1); % take first reachable node (not the souce node)


sort_perm = 1:N; % temp: assume tree is ordered. Otherwise: sort_perm = topological_sort(family_tree);

for i=1:N % loop on family members
    cur_member = sort_perm(i); cur_parents = parents(family_tree, cur_member);
    if(~isempty(cur_parents)) % only if node has parents
        K(cur_member, sort_perm(1:i-1)) = 0.5 .* ( ...
            K(cur_parents(1), sort_perm(1:i-1)) + K(cur_parents(2), sort_perm(1:i-1)) );
        K(sort_perm(1:i-1), cur_member) = K(cur_member, sort_perm(1:i-1))'; % make symmetric
        for j= [1 3] % loop on .. 
            U(cur_member, sort_perm(1:i-1),j) = ... % zero or two alleles shared
                ( U(cur_parents(1), sort_perm(1:i-1),j) + 0.5.*U(cur_parents(1), sort_perm(1:i-1),2) ) .* ...
                ( U(cur_parents(2), sort_perm(1:i-1),j) + 0.5.*U(cur_parents(2), sort_perm(1:i-1),2) );           
        end
    else
        U(cur_member, sort_perm(1:i-1), 1) = 1; U(cur_member, sort_perm(1:i-1), 3) = 0; % no allele sharing with new ancestors
    end
    U(cur_member, sort_perm(1:i-1), 2) = ...
        1 - U(cur_member, sort_perm(1:i-1), 1) - U(cur_member, sort_perm(1:i-1), 3); % one allele shared
    for j=1:3
        U(sort_perm(1:i-1), cur_member,j) = U(cur_member, sort_perm(1:i-1),j)'; % make symmetric
    end    
end % loop on family members
