% Test inference and mutual information computation for trees.

TOL = 0.00000000000000001;
n=5;  num_bnets = 1;

B = create_bounded_degree_bnet(n, 1, num_bnets);  



% Now transfer the bnets into my representation, and apply my way
BB = [];
for i=1:num_bnets
    BB{i} = B{i};
    % We need to switch the order ?
    for j=1:n
        BB{i}.rev_order(BB{i}.order(j)) = j;
    end
    

    % go over all the variables and give them their CPT's    
    B{i}.my_CPT = zeros(n,2);
    B{i}.my_CPT(B{i}.order(1),:) = CPD_to_CPT(B{i}.CPD{B{i}.order(1)})';
    for j=2:n
        tmp = CPD_to_CPT(B{i}.CPD{B{i}.order(j)})';
        B{i}.my_CPT(B{i}.order(j),:)  =  tmp(1,:);
    end

end


% First get the MI in the conventional way:
ttt = cputime; MI_mat = bnets_to_pairwise_MI(B,1);ttt_old_way = cputime - ttt
ttt = cputime; MY_MI_mat = bnets_to_pairwise_MI(B,1); ttt_my_way = cputime - ttt


% %     P = zuk_bnet_to_probs(B{i}); % Get the joint probability distribution
% %     % Compute all the marginals
% %     marginal_probs_mat = zeros(n,2);
% %     for j=1:n
% %         marginal_probs_mat(j,:) =  collapse_prob(P',j);
% %     end
% % 
% %     % Compute all the pairwise probabilities
% %     pairwise_probs_mat = zeros(n,n,4);
% %     for j=1:n
% %         for k=j+1:n
% %             pairwise_probs_mat(j,k,:) = collapse_prob(P',[j,k]);
% %         end
% %     end
% %     for j=1:4
% %         pairwise_probs_mat(:,:,j) = pairwise_probs_mat(:,:,j)+pairwise_probs_mat(:,:,j)';
% %     end
% % % SWAP 2 and 3 - will this help???
% %     tmp = pairwise_probs_mat(:,:,2);
% %     pairwise_probs_mat(:,:,2) = pairwise_probs_mat(:,:,3);
% %     pairwise_probs_mat(:,:,3) = tmp;
% % 
% %     T_Parents = zeros(1,n);
% %     for j=1:n
% %         tmp = find(BB{i}.dag(:,j));
% %         if(~isempty(tmp))
% %             T_Parents(j) = tmp;
% %         end
% %     end
% % pairwise_probs_mat = (pairwise_probs_mat + TOL) ./ (1 + 4*TOL);
% % MY_MI_mat2 = zeros(n,n);
% % MY_MI_mat2(:,:) = - repmat(marginal_probs_mat(:,1) .* log2(marginal_probs_mat(:,1)) + ...
% %     (1-marginal_probs_mat(:,1)) .* log2(1-marginal_probs_mat(:,1)), 1, n);
% % MY_MI_mat2 = MY_MI_mat2 + (1-1.*eye(n)).*MY_MI_mat2';
% % MY_MI_mat2 = MY_MI_mat2 + pairwise_probs_mat(:,:,1) .* log2(pairwise_probs_mat(:,:,1)) + ...
% %     pairwise_probs_mat(:,:,2) .* log2(pairwise_probs_mat(:,:,2)) + ...
% %     pairwise_probs_mat(:,:,3) .* log2(pairwise_probs_mat(:,:,3)) + ...
% %     pairwise_probs_mat(:,:,4) .* log2(pairwise_probs_mat(:,:,4));




max(max((MI_mat - MY_MI_mat).^2))




