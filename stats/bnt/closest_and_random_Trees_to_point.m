% Get the closest tree (in KL distance) to a point (distribution). 
% This can be done in polynomial time (up to the problem of representing P)
% 
% Input: 
% P - a probability distribution 
% 
% Output: 
% ClosestTree - tree closest to P 
% RandomTree - a randomly drawn tree
% MI_mat - mutual information matrix
% Closest_KL - Kullback-Leibler distance of closest tree to P 
% Random_KL - Kullback-Leibler distance of a random tree to P
% 
function [ClosestTree RandomTree MI_mat Closest_KL Random_KL] = ...
    closest_and_random_Trees_to_point(P)

n = log2(size(P,2));
num_probs = size(P,1);
% First compute mutual information matrix
MI_mat = zeros(n,n,num_probs);
% MI_mat2 = zeros(n,n,num_probs);

% first do a primitive calculation
for i=1:n
    for j=i+1:n
        MI_mat(i,j,:) = -entropy(collapse_prob(P,[i,j])');
   %%%     MI_mat2(i,j,:) = conditional_mutual_information(P, i,j,[]);
    end
end

for i=1:num_probs % heavy loop
    MI_mat(:,:,i) = MI_mat(:,:,i) + MI_mat(:,:,i)'; % make symmetric
end
marginal_entropy = zeros(n,num_probs);
for i=1:n
    marginal_entropy(i,:) = entropy(collapse_prob(P,i)');
end
MI_mat = MI_mat + ...
    reshape(repmat(reshape(marginal_entropy, 1, n, num_probs), 1, n), n, n, num_probs) + ...    
    reshape(repmat(reshape(marginal_entropy, 1, n, num_probs), n, 1), n, n, num_probs);




for i=1:num_probs % heavy loop
    MI_mat(:,:,i) = MI_mat(:,:,i) .* (1-eye(n)); %    MI_mat(:,:,i) = MI_mat(:,:,i) + MI_mat(:,:,i)'; % make mutual information symmetric

    ClosestTree = minimum_spanning_tree(-MI_mat(:,:,i)); % maximal spanning tree of MI matrix
    RandomTree = uniform_spanning_tree(n,0);
    
    Closest_KL = zeros(n,1); % Get relative entropy from the tree
    Closest_KL(i)=  -entropy(P(i,:)') -  sum(sum(MI_mat(:,:,i).*ClosestTree)) ./ 2;
    for j=1:n
        Closest_KL(i) = Closest_KL(i)+entropy(collapse_prob(P(i,:),j)');
    end
    Random_KL = zeros(n,1); 
    Random_KL(i)=  -entropy(P(i,:)') -  sum(sum(MI_mat(:,:,i).*RandomTree)) ./ 2;
    for j=1:n
        Random_KL(i) = Random_KL(i)+entropy(collapse_prob(P(i,:),j)');
    end

    
end

