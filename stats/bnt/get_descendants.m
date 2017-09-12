% Get all the descendant of X in G
% 
% Input: 
% X - set of nodes 
% G - directed acyclic graph 
% 
% Output: 
% des - set of X's descendants
% 
function des = get_descendants( X, G)

% Start with only X
des = X; old_des = [];
while (length(des) > length(old_des)) % As long as we reveal new vertices ...
    old_des = des;
    new_des =find(sum(G(des,:),1));
    des = union(des, new_des);
end
des = setdiff(des,X); % Remove the original vertex

return;