% Get the closest tree (in KL) to a point (distribution). 
% This can be done in polynomial time (up to the problem of representing P)
% 
% Input: 
% P - a probability distribution (might be more than one) 
%
% Output: 
% Tree 
% MI_mat 
% KL
% 
function [Tree MI_mat KL] = closest_Tree_to_point(P)

n = log2(size(P,2));
num_probs = size(P,1);

MI_mat = zeros(n,n,num_probs); % mutual information matrix
for i=1:n % first do a primitive calculation
    for j=i+1:n
        MI_mat(i,j,:) = conditional_mutual_information(P, i,j,[]);
    end
end
KL = zeros(num_probs,1); 
for i=1:num_probs % Here's the heavy loop
    MI_mat(:,:,i) = MI_mat(:,:,i) + MI_mat(:,:,i)'; % make symmetric
    Tree = minimum_spanning_tree(-MI_mat(:,:,i));  % find maximal spanning tree on the MI matrix    
    KL(i)=  -entropy(P(i,:)') -  sum(sum(MI_mat(:,:,i).*Tree)) ./ 2;  % Get relative entropy from tree
    for j=1:n
        KL(i) = KL(i)+entropy(collapse_prob(P(i,:),j)');
    end
end















