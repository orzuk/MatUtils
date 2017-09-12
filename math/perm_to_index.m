% Transferring a permutation to its lexicographical index
% 
% Input: 
% p - permutation on 1..N
% 
% Output: 
% ind - its index in all N! perumtations
% 
function ind = perm_to_index(p)
pp = p; num_perms = size(p,1); N = size(p,2);
ind = pp(:,1)*factorial(N-1);

for i=2:N
    pp = pp(:,2:end) - (pp(:,2:end) > repmat(pp(:,1), 1, N-i+1));
    ind = ind + pp(:,1)*factorial(N-i);
end

ind=ind+1;