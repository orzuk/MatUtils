% Get the closest point on a random tree (in KL) to a point (distribution). 
% This can be done in polynomial time (up to the problem of representing P)
function [Tree MI_mat KL] = random_Tree_to_point(P)

n = log2(size(P,2));
num_probs = size(P,1);
MI_mat = zeros(n,n,num_probs); % mutual information matrix
for i=1:n % first do a naive calculation
    for j=i+1:n
        MI_mat(i,j,:) = conditional_mutual_information(P, i,j,[]);
    end
end
KL = zeros(num_probs,1); 
for i=1:num_probs
    MI_mat(:,:,i) = MI_mat(:,:,i) + MI_mat(:,:,i)'; % make symmetric
    Tree = uniform_spanning_tree(n,0); % Draw a tree at 'random'
    KL(i)=  -entropy(P(i,:)') -  sum(sum(MI_mat(:,:,i).*Tree)) ./ 2;     % Get relative entropy from tree
    for j=1:n
        KL(i) = KL(i)+entropy(collapse_prob(P(i,:),j)');
    end
end

